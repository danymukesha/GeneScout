#' Calculate Codon Frequencies
#'
#' Calculate the frequency distribution of all 64 codons in a DNA sequence.
#'
#' @param sequence Character string or DNAString object representing DNA sequence
#' @param normalize Boolean indicating whether to normalize frequencies (default: TRUE)
#' @return Named numeric vector of codon frequencies
#' @export
#' @examples
#' sequence <- "ATGATGATGTTATTATTACGCCGCCGCC"
#' frequencies <- calculate_codon_frequencies(sequence)
calculate_codon_frequencies <- function(sequence, normalize = TRUE) {
    if (inherits(sequence, "DNAString")) {
        sequence <- Biostrings::toString(sequence)
    }

    if (!.validate_dna(sequence)) {
        stop("Invalid DNA sequence. Must contain only A, C, G, T (or N for unknown)")
    }

    adjusted_seq <- .adjust_codon_length(sequence, strategy = "truncate")

    if (nchar(adjusted_seq) < 3) {
        warning("Sequence too short for codon analysis (minimum 3 bp)")
        return(setNames(rep(0, 64), .get_all_codons()))
    }

    codons <- stringr::str_extract_all(adjusted_seq, ".{1,3}")[[1]]

    codon_table <- table(codons)
    all_codons <- .get_all_codons()

    full_table <- rep(0, 64)
    names(full_table) <- all_codons
    full_table[names(codon_table)] <- codon_table

    if (normalize) {
        total <- sum(full_table)
        if (total > 0) {
            full_table <- full_table / total
        }
    }

    return(full_table)
}

#' Calculate Shannon Entropy
#'
#' Calculate Shannon entropy of a codon frequency distribution.
#'
#' @param frequencies Named numeric vector of codon frequencies
#' @return Numeric value of Shannon entropy in bits
#' @export
#' @examples
#' sequence <- "ATGATGATGTTATTATTACGCCGCCGCC"
#' freqs <- calculate_codon_frequencies(sequence)
#' entropy <- calculate_shannon_entropy(freqs)
calculate_shannon_entropy <- function(frequencies) {
    frequencies <- frequencies[frequencies > 0]

    if (length(frequencies) == 0) {
        return(0)
    }

    entropy <- -sum(frequencies * log2(frequencies))
    return(entropy)
}

#' Calculate Kullback-Leibler Divergence
#'
#' Calculate the Kullback-Leibler (KL) divergence between two codon frequency
#' distributions. This measures how one distribution diverges from a second,
#' expected distribution.
#'
#' @param observed Named numeric vector of observed codon frequencies
#' @param reference Named numeric vector of reference codon frequencies
#' @param epsilon Small value to avoid log(0) (default: 1e-10)
#' @return Numeric value of KL divergence
#' @export
#' @examples
#' sequence1 <- "ATGATGATGTTATTATTACGCCGCCGCC"
#' sequence2 <- "TTATTATTACGCCGCCGCCATGATGATG"
#' freqs1 <- calculate_codon_frequencies(sequence1)
#' freqs2 <- calculate_codon_frequencies(sequence2)
#' kl_div <- calculate_kl_divergence(freqs1, freqs2)
calculate_kl_divergence <- function(observed, reference, epsilon = 1e-10) {
    if (!all(names(observed) == names(reference))) {
        stop("Observed and reference frequencies must have the same codon order")
    }

    kl_sum <- 0
    for (i in seq_along(observed)) {
        p <- max(observed[i], epsilon)
        q <- max(reference[i], epsilon)

        if (p > 0 && q > 0) {
            kl_sum <- kl_sum + p * log2(p / q)
        }
    }

    return(kl_sum)
}

#' Calculate Relative Synonymous Codon Usage (RSCU)
#'
#' Calculate the Relative Synonymous Codon Usage for a codon frequency distribution.
#' RSCU > 1 indicates codon usage bias towards that codon.
#'
#' @param frequencies Named numeric vector of codon frequencies
#' @return Named numeric vector of RSCU values
#' @export
#' @examples
#' sequence <- "ATGATGATGTTATTATTACGCCGCCGCC"
#' freqs <- calculate_codon_frequencies(sequence)
#' rscu <- calculate_rscu(freqs)
calculate_rscu <- function(frequencies) {
    all_codons <- .get_all_codons()

    aa_table <- c(
        "ATA" = "I", "ATC" = "I", "ATT" = "I", "ATG" = "M",
        "ACA" = "T", "ACC" = "T", "ACG" = "T", "ACT" = "T",
        "AAC" = "N", "AAT" = "N", "AAA" = "K", "AAG" = "K",
        "AGC" = "S", "AGT" = "S", "AGA" = "R", "AGG" = "R",
        "CTA" = "L", "CTC" = "L", "CTG" = "L", "CTT" = "L",
        "CCA" = "P", "CCC" = "P", "CCG" = "P", "CCT" = "P",
        "CAC" = "H", "CAT" = "H", "CAA" = "Q", "CAG" = "Q",
        "CGA" = "R", "CGC" = "R", "CGG" = "R", "CGT" = "R",
        "GTA" = "V", "GTC" = "V", "GTG" = "V", "GTT" = "V",
        "GCA" = "A", "GCC" = "A", "GCG" = "A", "GCT" = "A",
        "GAC" = "D", "GAT" = "D", "GAA" = "E", "GAG" = "E",
        "GGA" = "G", "GGC" = "G", "GGG" = "G", "GGT" = "G",
        "TCA" = "S", "TCC" = "S", "TCG" = "S", "TCT" = "S",
        "TTC" = "F", "TTT" = "F", "TTA" = "L", "TTG" = "L",
        "TAC" = "Y", "TAT" = "Y", "TAA" = "*", "TAG" = "*",
        "TGC" = "C", "TGT" = "C", "TGA" = "*", "TGG" = "W"
    )

    aa_names <- aa_table[all_codons]
    unique_aa <- unique(aa_names)

    rscu_values <- numeric(64)
    names(rscu_values) <- all_codons

    for (aa in unique_aa) {
        if (aa == "*") next

        aa_codons <- names(aa_names[aa_names == aa])
        aa_freqs <- frequencies[aa_codons]
        aa_sum <- sum(aa_freqs)

        if (aa_sum > 0) {
            n_codons <- length(aa_codons)
            expected_freq <- aa_sum / n_codons
            rscu_values[aa_codons] <- aa_freqs / expected_freq
        }
    }

    return(rscu_values)
}

#' Calculate Effective Number of Codons (ENC)
#'
#' Calculate the Effective Number of Codons (ENC), a measure of codon usage bias.
#' ENC ranges from 20 (extreme bias) to 61 (no bias).
#'
#' @param frequencies Named numeric vector of codon frequencies
#' @return Numeric value of ENC
#' @export
#' @examples
#' sequence <- "ATGATGATGTTATTATTACGCCGCCGCC"
#' freqs <- calculate_codon_frequencies(sequence)
#' enc <- calculate_enc(freqs)
calculate_enc <- function(frequencies) {
    aa_table <- c(
        "ATA" = "I", "ATC" = "I", "ATT" = "I",
        "ACA" = "T", "ACC" = "T", "ACG" = "T", "ACT" = "T",
        "AAC" = "N", "AAT" = "N", "AAA" = "K", "AAG" = "K",
        "AGC" = "S", "AGT" = "S", "AGA" = "R", "AGG" = "R",
        "CTA" = "L", "CTC" = "L", "CTG" = "L", "CTT" = "L",
        "CCA" = "P", "CCC" = "P", "CCG" = "P", "CCT" = "P",
        "CAC" = "H", "CAT" = "H", "CAA" = "Q", "CAG" = "Q",
        "CGA" = "R", "CGC" = "R", "CGG" = "R", "CGT" = "R",
        "GTA" = "V", "GTC" = "V", "GTG" = "V", "GTT" = "V",
        "GCA" = "A", "GCC" = "A", "GCG" = "A", "GCT" = "A",
        "GAC" = "D", "GAT" = "D", "GAA" = "E", "GAG" = "E",
        "GGA" = "G", "GGC" = "G", "GGG" = "G", "GGT" = "G",
        "TCA" = "S", "TCC" = "S", "TCG" = "S", "TCT" = "S",
        "TTC" = "F", "TTT" = "F", "TTA" = "L", "TTG" = "L",
        "TAC" = "Y", "TAT" = "Y",
        "TGC" = "C", "TGT" = "C",
        "TGG" = "W", "ATG" = "M"
    )

    all_codons <- .get_all_codons()
    aa_names <- aa_table[all_codons]
    unique_aa <- unique(aa_names)

    enc_sum <- 0
    n_aa <- 0

    for (aa in unique_aa) {
        if (aa == "*") next

        aa_codons <- names(aa_names[aa_names == aa])
        aa_freqs <- frequencies[aa_codons]
        aa_freqs <- aa_freqs[aa_freqs > 0]

        if (length(aa_freqs) > 0) {
            aa_entropy <- -sum(aa_freqs / sum(aa_freqs) * log2(aa_freqs / sum(aa_freqs)))
            n_codons <- length(aa_codons)

            if (n_codons > 1) {
                enc_sum <- enc_sum + aa_entropy / log2(n_codons)
                n_aa <- n_aa + 1
            }
        }
    }

    if (n_aa > 0) {
        enc <- 2 + 9 / enc_sum
        return(min(enc, 61))
    } else {
        return(61)
    }
}
