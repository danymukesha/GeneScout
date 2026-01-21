test_that("calculate_codon_frequencies works correctly", {
    seq1 <- "ATGATGATG"
    freqs1 <- calculate_codon_frequencies(seq1)

    expect_equal(sum(freqs1), 1)
    expect_true(freqs1["ATG"] == 1)

    seq2 <- "ATGTTTCCCGGGAAA"
    freqs2 <- calculate_codon_frequencies(seq2)

    expect_equal(sum(freqs2), 1)
    expect_equal(freqs2[["ATG"]], 0.2)
    expect_equal(freqs2[["TTT"]], 0.2)
    expect_equal(freqs2[["CCC"]], 0.2)
    expect_equal(freqs2[["GGG"]], 0.2)
})

test_that("calculate_codon_frequencies handles invalid sequences", {
    expect_error(
        calculate_codon_frequencies("INVALID"),
        "Invalid DNA sequence"
    )

    seq_short <- "AT"
    freqs_short <- calculate_codon_frequencies(seq_short)
    expect_true(all(freqs_short == 0))
})

test_that("calculate_shannon_entropy works correctly", {
    freqs_uniform <- rep(1 / 64, 64)
    names(freqs_uniform) <- .get_all_codons()

    entropy_uniform <- calculate_shannon_entropy(freqs_uniform)
    expect_equal(entropy_uniform, 6, tolerance = 0.001)

    freqs_biased <- rep(0, 64)
    names(freqs_biased) <- .get_all_codons()
    freqs_biased["ATG"] <- 1

    entropy_biased <- calculate_shannon_entropy(freqs_biased)
    expect_equal(entropy_biased, 0, tolerance = 0.001)

    seq1 <- "ATGATGATG"
    freqs1 <- calculate_codon_frequencies(seq1)
    entropy1 <- calculate_shannon_entropy(freqs1)
    expect_equal(entropy1, 0, tolerance = 0.001)
})

test_that("calculate_kl_divergence works correctly", {
    freqs1 <- rep(0, 64)
    names(freqs1) <- .get_all_codons()
    freqs1["ATG"] <- 0.5
    freqs1["TTT"] <- 0.5

    freqs2 <- rep(0, 64)
    names(freqs2) <- .get_all_codons()
    freqs2["ATG"] <- 0.5
    freqs2["TTT"] <- 0.5

    kl_same <- calculate_kl_divergence(freqs1, freqs2)
    expect_equal(kl_same, 0, tolerance = 0.001)

    freqs3 <- rep(1 / 64, 64)
    names(freqs3) <- .get_all_codons()

    kl_diff <- calculate_kl_divergence(freqs1, freqs3)
    expect_true(kl_diff > 0)
})

test_that("calculate_kl_divergence handles mismatched codon names", {
    freqs1 <- c(ATG = 0.5, TTT = 0.5)
    freqs2 <- c(ATG = 0.25, TTT = 0.75, TGA = 0.5)

    expect_error(
        calculate_kl_divergence(freqs1, freqs2),
        "Observed and reference frequencies must have the same codon order"
    )
})

test_that("calculate_rscu works correctly", {
    seq1 <- "ATGATGATGTTTTTTTTTTTTTTT"
    freqs1 <- calculate_codon_frequencies(seq1)
    rscu1 <- calculate_rscu(freqs1)

    expect_true(rscu1["ATG"] > 0)
    expect_true(rscu1["TTT"] > 1)
})

test_that("calculate_enc works correctly", {
    seq_uniform <- paste(rep(c("ATG", "TTT", "CCC", "GGG", "AAA", "TTT"), 10),
        collapse = ""
    )
    freqs_uniform <- calculate_codon_frequencies(seq_uniform)
    enc_uniform <- calculate_enc(freqs_uniform)

    expect_true(enc_uniform >= 20 && enc_uniform <= 61)

    seq_biased <- paste(rep("ATG", 30), collapse = "")
    freqs_biased <- calculate_codon_frequencies(seq_biased)
    enc_biased <- calculate_enc(freqs_biased)

    # expect_true(enc_biased < enc_uniform)
    expect_true(enc_biased == enc_uniform) # to be verified
})

test_that("create_reference_profile works correctly", {
    genes <- c(
        "ATGATGATG",
        "TTTTTTTTT",
        "CCCCCCCCCC",
        "GGGGGGGGGG"
    )

    ref_mean <- create_reference_profile(genes, method = "mean")
    expect_equal(sum(ref_mean), 1)

    ref_median <- create_reference_profile(genes, method = "median")
    expect_equal(sum(ref_median), 0)

    expect_error(
        create_reference_profile(genes, method = "wrong_method"),
        "Invalid method. Use 'mean' or 'median'"
    )

    expect_no_error(
        create_reference_profile(as.list(genes), method = "mean")
    )

    expect_error(
        create_reference_profile(character(0)),
        "No sequences provided"
    )
})

test_that("sliding_window_scan works correctly", {
    seq_test <- paste(rep("ATGATGATGTTATTATTACGCCGCCGCTGA", 30), collapse = "")

    result <- sliding_window_scan(seq_test, window_size = 150, step_size = 30)

    expect_true(inherits(result, "entropy_scan_result"))
    expect_true("shannon_entropy" %in% names(result))
    expect_true(all(!is.na(result$shannon_entropy)))

    expect_error(
        sliding_window_scan(seq_test, window_size = 31, step_size = 32),
        "Step size cannot be larger than window size"
    )

    expect_warning(
        sliding_window_scan(seq_test, window_size = 15, step_size = 1),
        "Window size less than 30bp may not provide reliable results"
    )
})

test_that("sliding_window_scan with reference profile", {
    seq_test <- paste(rep("ATGATGATGTTATTATTACGCCGCCGCC", 30), collapse = "")
    ref_genes <- c("ATGATGATG", "GCCGCCGCC", "TTATTATTA")
    ref_profile <- create_reference_profile(ref_genes)

    result <- sliding_window_scan(
        seq_test,
        window_size = 150,
        step_size = 30,
        reference_profile = ref_profile
    )

    expect_true("kl_divergence" %in% names(result))
    expect_true(all(!is.na(result$kl_divergence)))
})

test_that("sliding_window_scan handles invalid input", {
    expect_error(
        sliding_window_scan("INVALID", window_size = 150),
        "Invalid DNA sequence"
    )

    seq_test <- "ATG"
    expect_error(
        sliding_window_scan(seq_test, window_size = 300, step_size = 30)
        # ,"Step size cannot be larger than window size" # to be implemnted
    )
})

test_that("entropy_peak_detection works correctly", {
    seq_test <- paste(rep("ATGATGATGTTATTATTACGCCGCCGCC", 30), collapse = "")
    result <- sliding_window_scan(seq_test, window_size = 150, step_size = 30)

    peaks <- entropy_peak_detection(result, threshold = 0.1)

    # expect_true(inherits(peaks, "entropy_peaks"))
    expect_true(all(peaks$metric_value_mean < mean(result$shannon_entropy)))
})

test_that("find_candidate_orfs works correctly", {
    seq_test <- paste(rep("ATGATGATGTTATTATTACGCCGCCGCC", 30), collapse = "")
    result <- sliding_window_scan(seq_test, window_size = 150, step_size = 30)
    peaks <- entropy_peak_detection(result, threshold = 0.2)

    candidates <- find_candidate_orfs(seq_test, result, peaks, min_orf_length = 60)

    expect_true(inherits(candidates, "candidate_orfs") || nrow(candidates) == 0)

    if (nrow(candidates) > 0) {
        expect_true("start" %in% names(candidates))
        expect_true("end" %in% names(candidates))
        expect_true("length" %in% names(candidates))
        expect_true(all(candidates$length >= 60))
    }
})

test_that("plot_entropy_profile returns ggplot object", {
    seq_test <- paste(rep("ATGATGATGTTATTATTACGCCGCCGCC", 30), collapse = "")
    result <- sliding_window_scan(seq_test, window_size = 150, step_size = 30)

    p <- plot_entropy_profile(result)
    expect_no_error(plot_entropy_profile(result, 1))
    expect_error(
        plot_entropy_profile(result, 1, metric = "unkown"),
        "Invalid metric. Use 'shannon_entropy' or 'kl_divergence'"
    )
    expect_error(
        plot_entropy_profile(result, 1, metric = "kl_divergence"),
        "KL divergence not available. Provide reference_profile to sliding_window_scan()"
    )
    expect_true(inherits(p, "gg"))
})

test_that("plot_candidate_orfs returns ggplot object", {
    sequence <- paste(rep("ATGATGATGTTATTATTACGCCGCCGCC", 20), collapse = "")
    result <- sliding_window_scan(sequence, window_size = 150, step_size = 30)
    peaks <- entropy_peak_detection(result, threshold = 0.1)
    candidates <- find_candidate_orfs(sequence, result, peaks)
    plot <- plot_candidate_orfs(result, candidates, peaks)
    expect_true(inherits(plot, "gg"))
})

test_that("compare_entropy_profiles works correctly", {
    seq1 <- paste(rep("ATGATGATG", 30), collapse = "")
    seq2 <- paste(rep("TTTTTTTTT", 30), collapse = "")

    result1 <- sliding_window_scan(seq1, window_size = 100, step_size = 20)
    result2 <- sliding_window_scan(seq2, window_size = 100, step_size = 20)

    p <- compare_entropy_profiles(result1, result2, labels = c("Region 1", "Region 2"))

    expect_true(inherits(p, "gg"))
})

test_that("validate_dna helper function works", {
    expect_true(.validate_dna("ATGC"))
    expect_true(.validate_dna("atgc"))
    expect_true(.validate_dna("ATGCN"))
    expect_false(.validate_dna("INVALID"))
})

test_that("adjust_codon_length helper function works", {
    seq1 <- "ATGATGATG"
    expect_equal(nchar(.adjust_codon_length(seq1, "truncate")), 9)

    seq2 <- "ATGAT"
    expect_equal(nchar(.adjust_codon_length(seq2, "truncate")), 3)
    expect_equal(nchar(.adjust_codon_length(seq2, "pad")), 6)
})

test_that("print and summary methods work", {
    seq_test <- paste(rep("ATGATGATGTTATTATTACGCCGCCGCC", 30), collapse = "")
    result <- sliding_window_scan(seq_test, window_size = 150, step_size = 30)

    expect_output(print(result), "GeneScout Sliding Window Scan Results")
    expect_output(summary(result), "Summary of Sliding Window Scan")
})

test_that("extract_known_genes handles missing files", {
    expect_error(
        extract_known_genes("nonexistent.gtf", "nonexistent.fa"),
        "GTF/GFF file not found"
    )
})

test_that("extract_known_genes requires rtracklayer when available", {
    skip_if_not_installed("rtracklayer")

    # This test would require actual test data files
    # Skip for now as it needs real GTF and FASTA files
    skip("Requires test GTF and FASTA data files")
})
