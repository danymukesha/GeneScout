#' Create Reference Codon Profile
#'
#' Create a reference codon usage profile from known genes or a set of training sequences.
#' This profile represents the organism's codon usage bias pattern.
#'
#' @param sequences Character vector, DNAStringSet, or list of DNA sequences
#' @param method Character string: "mean" or "median" for aggregating frequencies
#' @return Named numeric vector representing reference codon frequencies
#' @export
#' @examples
#' known_genes <- c("ATGATGATG", "GCCGCCGCC", "TTATTATTA")
#' ref_profile <- create_reference_profile(known_genes)
create_reference_profile <- function(sequences, method = "mean") {
    if (inherits(sequences, "DNAStringSet")) {
        sequences <- as.character(sequences)
    } else if (is.list(sequences)) {
        sequences <- unlist(sequences)
    }

    if (length(sequences) == 0) {
        stop("No sequences provided")
    }

    freq_matrix <- sapply(sequences, function(seq) {
        calculate_codon_frequencies(seq, normalize = TRUE)
    })

    if (method == "mean") {
        ref_profile <- rowMeans(freq_matrix)
    } else if (method == "median") {
        ref_profile <- apply(freq_matrix, 1, median)
    } else {
        stop("Invalid method. Use 'mean' or 'median'")
    }

    ref_profile <- ref_profile / sum(ref_profile)
    return(ref_profile)
}

#' Sliding Window Entropy Scan
#'
#' Perform a sliding window scan across a DNA sequence to calculate entropy and
#' KL divergence metrics at each position. This is the main function for de novo
#' ORF discovery.
#'
#' @param sequence Character string or DNAString object representing DNA sequence
#' @param window_size Integer size of sliding window in base pairs (default: 300)
#' @param step_size Integer step size for sliding window (default: 30)
#' @param reference_profile Optional reference codon frequency profile. If NULL,
#'   only Shannon entropy is calculated.
#' @param min_codons Minimum number of complete codons required for analysis (default: 10)
#' @return A data frame with window coordinates and entropy metrics
#' @export
#' @examples
#' sequence <- paste(rep("ATGATGATGTTATTATTACGCCGCCGCC", 20), collapse = "")
#' result <- sliding_window_scan(sequence, window_size = 150, step_size = 30)
sliding_window_scan <- function(sequence,
                                window_size = 300,
                                step_size = 30,
                                reference_profile = NULL,
                                min_codons = 10) {
    if (inherits(sequence, "DNAString")) {
        sequence <- Biostrings::toString(sequence)
    }

    if (!.validate_dna(sequence)) {
        stop("Invalid DNA sequence. Must contain only A, C, G, T (or N for unknown)")
    }

    seq_length <- nchar(sequence)

    if (window_size < 30) {
        warning("Window size less than 30bp may not provide reliable results")
    }

    if (step_size > window_size) {
        stop("Step size cannot be larger than window size")
    }

    window_coords <- .get_window_coords(seq_length, window_size, step_size)

    # Store all incidences for final report
    mismatch_report <- list()

    results <- lapply(1:nrow(window_coords), function(i) {
        start <- window_coords$start[i]
        end <- window_coords$end[i]
        window_seq <- substr(sequence, start, end)

        codon_freqs <- calculate_codon_frequencies(window_seq, normalize = TRUE)
        shannon <- calculate_shannon_entropy(codon_freqs)

        result_row <- data.frame(
            window_id = window_coords$window_id[i],
            start = start,
            end = end,
            shannon_entropy = shannon,
            kl_divergence = NA_real_,
            stringsAsFactors = FALSE
        )

        if (!is.null(reference_profile)) {
            extra_codons <- setdiff(names(codon_freqs), names(reference_profile))
            missing_codons <- setdiff(names(reference_profile), names(codon_freqs))

            if (length(extra_codons) > 0 || length(missing_codons) > 0) {
                # Save mismatch info instead of printing now
                mismatch_report[[length(mismatch_report) + 1]] <<- list(
                    window_id = window_coords$window_id[i],
                    extra_codons = extra_codons,
                    missing_codons = missing_codons
                )
            }

            # Align codons and fill missing ones
            aligned_codons <- intersect(names(reference_profile), names(codon_freqs))
            codon_freqs_aligned <- codon_freqs[aligned_codons]
            reference_profile_aligned <- reference_profile[aligned_codons]

            for (codon in setdiff(names(reference_profile), names(codon_freqs_aligned))) {
                codon_freqs_aligned[codon] <- 0
            }

            codon_freqs_aligned <- codon_freqs_aligned[names(reference_profile)]
            result_row$kl_divergence <- calculate_kl_divergence(codon_freqs_aligned, reference_profile)
        }

        return(result_row)
    })

    results_df <- dplyr::bind_rows(results)
    class(results_df) <- c("entropy_scan_result", class(results_df))

    # FINAL REPORT for codon mismatches
    if (length(mismatch_report) > 0) {
        message("⚠ Codon mismatch summary for this sequence:")
        for (entry in mismatch_report) {
            msg <- paste0(
                "  Window ", entry$window_id, ": ",
                if (length(entry$extra_codons) > 0) paste0("Extra codons = ", paste(entry$extra_codons, collapse = ", ")) else "",
                if (length(entry$missing_codons) > 0) {
                    paste0(
                        if (length(entry$extra_codons) > 0) "; " else "",
                        "Missing codons = ", paste(entry$missing_codons, collapse = ", ")
                    )
                } else {
                    ""
                }
            )
            message(msg)
        }
    }

    return(results_df)
}


#' @export
print.entropy_scan_result <- function(x, ...) {
    cat("GeneScout Sliding Window Scan Results\n")
    cat("=======================================\n")
    cat("Number of windows:", nrow(x), "\n")
    cat("Sequence range:", x$start[1], "-", x$end[nrow(x)], "bp\n")
    cat("\nEntropy Statistics:\n")
    cat("  Mean Shannon Entropy:", round(mean(x$shannon_entropy), 3), "bits\n")
    cat("  Std. Dev. Shannon Entropy:", round(sd(x$shannon_entropy), 3), "bits\n")
    cat("  Min Shannon Entropy:", round(min(x$shannon_entropy), 3), "bits\n")
    cat("  Max Shannon Entropy:", round(max(x$shannon_entropy), 3), "bits\n")

    if (!is.na(x$kl_divergence[1])) {
        cat("\nKL Divergence Statistics:\n")
        cat("  Mean KL Divergence:", round(mean(x$kl_divergence, na.rm = TRUE), 3), "\n")
        cat("  Std. Dev. KL Divergence:", round(sd(x$kl_divergence, na.rm = TRUE), 3), "\n")
    }

    invisible(x)
}

#' @export
summary.entropy_scan_result <- function(object, ...) {
    cat("Summary of Sliding Window Scan\n")
    cat("==============================\n\n")

    entropy_stats <- summary(object$shannon_entropy)
    print(entropy_stats)

    if (!all(is.na(object$kl_divergence))) {
        cat("\nKL Divergence:\n")
        kl_stats <- summary(object$kl_divergence, na.rm = TRUE)
        print(kl_stats)
    }

    invisible(object)
}

# Future implementations:
##  1. Performance Optimization (The "C++" Factor)
### Scanning a whole chromosome (millions of letters) with a sliding window
### in pure R will be slow.
### - The fix: Considering using the `Rcpp` package to move
### the `sliding_window_scan` loop into C++. Since you are just doing simple
###  math on strings, C++ will make it **100x faster** (just an estimation).

## 2. Reverse Complement Scanning
### Currently, your code only scans the "Forward" strand of DNA. In biology,
### genes can be on the "Reverse" strand (the DNA runs the other way).
### - Update `sliding_window_scan` to also scan the reverse complement.
# ```r
# # Simple logic to add:
# rev_seq <- Biostrings::reverseComplement(Biostrings::DNAString(sequence))
# # Run your scan on both and label them "Forward" and "Reverse"
# ```

##  3. Better Peak Detection
### The `entropy_peak_detection` uses a simple quantile threshold.
### In genomics, "noise" changes across the chromosome.
### - The fix: Use a **Z-score** or a **Moving Average** to find peaks.
### This helps find genes in "high-GC" areas where the baseline entropy
### might be naturally different.

##  4. Handling "N" (Unknown Bases)
### Add the validator checks for `N`, but if a window has 50% `N`s,
### the entropy calculation will be misleading.
### - The fix: Add a parameter `max_n_threshold = 0.1`. If a window has more
### than 10% "N" bases, skip it or mark it as "Low Confidence."
