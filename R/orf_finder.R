#' Detect Entropy Peaks
#'
#' Identify peaks in entropy profiles that may indicate potential coding regions.
#'
#' @param scan_result Data frame from sliding_window_scan()
#' @param metric Character string: "shannon_entropy" or "kl_divergence"
#' @param method Character string: "quantile" or "sd" for threshold calculation
#' @param threshold Quantile threshold (for "quantile" method) or number of SDs (for "sd" method)
#' @param min_peak_width Minimum width of peaks in windows (default: 3)
#' @return Data frame with peak information
#' @export
#' @examples
#' sequence <- paste(rep("ATGATGATGTTATTATTACGCCGCCGCC", 20), collapse = "")
#' result <- sliding_window_scan(sequence, window_size = 150, step_size = 30)
#' peaks <- entropy_peak_detection(result, threshold = 0.1)
entropy_peak_detection <- function(scan_result,
                                   metric = "shannon_entropy",
                                   method = "quantile",
                                   threshold = 0.1,
                                   min_peak_width = 3) {
    if (!inherits(scan_result, "entropy_scan_result")) {
        scan_result <- structure(scan_result, class = c("entropy_scan_result", class(scan_result)))
    }

    if (!metric %in% c("shannon_entropy", "kl_divergence")) {
        stop("Invalid metric. Use 'shannon_entropy' or 'kl_divergence'")
    }

    if (metric == "kl_divergence" && all(is.na(scan_result$kl_divergence))) {
        stop("KL divergence not available. Provide reference_profile to sliding_window_scan()")
    }

    values <- scan_result[[metric]]

    if (method == "quantile") {
        cutoff <- quantile(values, threshold, na.rm = TRUE)
    } else if (method == "sd") {
        cutoff <- mean(values, na.rm = TRUE) - threshold * sd(values, na.rm = TRUE)
        if (metric == "kl_divergence") {
            cutoff <- mean(values, na.rm = TRUE) + threshold * sd(values, na.rm = TRUE)
        }
    } else {
        stop("Invalid method. Use 'quantile' or 'sd'")
    }

    below_threshold <- values <= cutoff
    if (metric == "kl_divergence") {
        below_threshold <- values >= cutoff
    }

    peak_groups <- rle(below_threshold)
    peak_starts <- cumsum(c(1, peak_groups$lengths[-length(peak_groups$lengths)]))
    peak_starts <- peak_starts[peak_groups$values]
    peak_lengths <- peak_groups$lengths[peak_groups$values]

    valid_peaks <- peak_lengths >= min_peak_width
    peak_starts <- peak_starts[valid_peaks]
    peak_lengths <- peak_lengths[valid_peaks]

    if (length(peak_starts) == 0) {
        warning("No peaks detected with current threshold")
        return(data.frame())
    }

    peaks_df <- data.frame(
        peak_id = 1:length(peak_starts),
        start_window = peak_starts,
        end_window = peak_starts + peak_lengths - 1,
        num_windows = peak_lengths,
        stringsAsFactors = FALSE
    )

    peaks_df <- dplyr::mutate(peaks_df,
        start_bp = scan_result$start[start_window],
        end_bp = scan_result$end[end_window],
        metric_value_mean = sapply(1:nrow(peaks_df), function(i) {
            mean(scan_result[[metric]][peaks_df$start_window[i]:peaks_df$end_window[i]], na.rm = TRUE)
        }),
        metric_value_min = sapply(1:nrow(peaks_df), function(i) {
            min(scan_result[[metric]][peaks_df$start_window[i]:peaks_df$end_window[i]], na.rm = TRUE)
        })
    )

    class(peaks_df) <- c("entropy_peaks", class(peaks_df))
    return(peaks_df)
}

#' Find Candidate ORFs
#'
#' Identify candidate open reading frames (ORFs) based on entropy analysis and
#' sequence features.
#'
#' @param sequence Character string or DNAString object representing DNA sequence
#' @param scan_result Data frame from sliding_window_scan()
#' @param peaks Data frame from entropy_peak_detection()
#' @param min_orf_length Minimum ORF length in base pairs (default: 90, i.e., 30 aa)
#' @param start_codons Character vector of valid start codons (default: c("ATG"))
#' @param stop_codons Character vector of stop codons (default: c("TAA", "TAG", "TGA"))
#' @return Data frame with candidate ORF information
#' @export
#' @examples
#' sequence <- paste(rep("ATGATGATGTTATTATTACGCCGCCGCC", 20), collapse = "")
#' result <- sliding_window_scan(sequence, window_size = 150, step_size = 30)
#' peaks <- entropy_peak_detection(result, threshold = 0.1)
#' candidates <- find_candidate_orfs(sequence, result, peaks)
find_candidate_orfs <- function(sequence,
                                scan_result,
                                peaks,
                                min_orf_length = 90,
                                start_codons = c("ATG"),
                                stop_codons = c("TAA", "TAG", "TGA")) {
    if (inherits(sequence, "DNAString")) {
        sequence <- Biostrings::toString(sequence)
    }

    if (nrow(peaks) == 0) {
        warning("No peaks provided, returning empty data frame")
        return(data.frame())
    }

    candidates <- lapply(1:nrow(peaks), function(i) {
        peak_start <- peaks$start_bp[i]
        peak_end <- peaks$end_bp[i]

        peak_seq <- substr(sequence, peak_start, peak_end)

        found_orfs <- find_orfs_in_sequence(
            peak_seq,
            min_orf_length = min_orf_length,
            start_codons = start_codons,
            stop_codons = stop_codons
        )

        if (nrow(found_orfs) == 0) {
            return(NULL)
        }

        found_orfs$start <- found_orfs$start + peak_start - 1
        found_orfs$end <- found_orfs$end + peak_start - 1
        found_orfs$peak_id <- i
        found_orfs$entropy_score <- peaks$metric_value_min[i]

        return(found_orfs)
    })

    candidates_df <- dplyr::bind_rows(candidates)

    if (nrow(candidates_df) > 0) {
        candidates_df <- dplyr::arrange(candidates_df, start)
        class(candidates_df) <- c("candidate_orfs", class(candidates_df))
    }

    return(candidates_df)
}

#' Find ORFs in Sequence
#'
#' Helper function to find ORFs within a specific sequence region.
#'
#' @param sequence Character string of DNA sequence
#' @param min_orf_length Minimum ORF length in base pairs
#' @param start_codons Character vector of start codons
#' @param stop_codons Character vector of stop codons
#' @return Data frame with ORF information
#' @keywords internal
find_orfs_in_sequence <- function(sequence,
                                  min_orf_length,
                                  start_codons,
                                  stop_codons) {
    seq_len <- nchar(sequence)

    start_positions <- c()

    for (codon in start_codons) {
        matches <- gregexpr(codon, sequence)[[1]]
        if (matches[1] != -1) {
            start_positions <- c(start_positions, matches)
        }
    }

    start_positions <- unique(start_positions)
    start_positions <- sort(start_positions)

    orfs <- lapply(start_positions, function(start_pos) {
        seq_from_start <- substr(sequence, start_pos, seq_len)

        if (nchar(seq_from_start) < min_orf_length) {
            return(NULL)
        }

        stop_pos <- NA

        for (stop_codon in stop_codons) {
            matches <- gregexpr(stop_codon, seq_from_start)[[1]]
            if (matches[1] != -1) {
                valid_stops <- matches[matches %% 3 == 0]
                if (length(valid_stops) > 0) {
                    candidate <- start_pos + valid_stops[1] + 2
                    if (candidate - start_pos + 1 >= min_orf_length) {
                        if (is.na(stop_pos) || candidate < stop_pos) {
                            stop_pos <- candidate
                        }
                    }
                }
            }
        }

        if (is.na(stop_pos)) {
            orf_seq <- seq_from_start
            stop_codon_found <- FALSE
            orf_end <- start_pos + nchar(orf_seq) - 1
        } else {
            orf_seq <- substr(sequence, start_pos, stop_pos)
            stop_codon_found <- TRUE
            orf_end <- stop_pos
        }

        if (nchar(orf_seq) >= min_orf_length) {
            return(data.frame(
                start = start_pos,
                end = orf_end,
                length = nchar(orf_seq),
                sequence = orf_seq,
                has_stop_codon = stop_codon_found,
                frame = ((start_pos - 1) %% 3) + 1,
                stringsAsFactors = FALSE
            ))
        } else {
            return(NULL)
        }
    })

    orfs_df <- dplyr::bind_rows(orfs)

    if (!is.null(orfs_df) && nrow(orfs_df) > 0) {
        orfs_df <- orfs_df[!duplicated(orfs_df$start), ]
        orfs_df <- dplyr::arrange(orfs_df, start)
    }

    return(orfs_df)
}
