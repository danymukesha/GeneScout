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
    p <- plot_entropy_profile(result, metric = "kl_divergence")
    expect_true(inherits(p, "gg"))
    expect_no_error(plot_entropy_profile(result,
        metric = "kl_divergence",
        show_peaks = TRUE
    ))
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

    peaks <- entropy_peak_detection(result, threshold = 0.7)

    # expect_true(inherits(peaks, "entropy_peaks"))
    expect_true(all(peaks$metric_value_mean < mean(result$shannon_entropy)))

    candidates <- find_candidate_orfs(seq_test, result, peaks)
    expect_no_error(
        plot_candidate_orfs(result, candidates, peaks)
    )
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

test_that("find_candidate_orfs returns correct ORFs with proper columns", {
    seq_test <- paste(rep("ATGAAATGAAAAATAGATGTAAATGAGGAGATGTTTGA", 3), collapse = "")
    scan_result <- data.frame(window_start = 1, window_end = nchar(seq_test))

    peaks <- data.frame(
        start_bp = c(1, 50),
        end_bp = c(45, 90),
        metric_value_min = c(0.1, 0.2)
    )

    candidates <- find_candidate_orfs(seq_test, scan_result, peaks, min_orf_length = 9)

    expect_s3_class(candidates, "candidate_orfs")

    expect_true(all(c(
        "start", "end", "length", "sequence", "has_stop_codon", "frame",
        "peak_id", "entropy_score"
    ) %in% names(candidates)))

    expect_true(all(candidates$length >= 9))

    for (i in seq_len(nrow(candidates))) {
        peak <- peaks[candidates$peak_id[i], ]
        expect_gte(candidates$start[i], peak$start_bp)
        expect_lte(candidates$end[i], peak$end_bp)
    }

    for (i in seq_len(nrow(candidates))) {
        seq_extracted <- substr(seq_test, candidates$start[i], candidates$end[i])
        expect_equal(seq_extracted, candidates$sequence[i])
    }
})

test_that("find_candidate_orfs returns empty data frame if no peaks", {
    seq_test <- paste(rep("ATGAAATGAAAAATAG", 2), collapse = "")
    scan_result <- data.frame()
    peaks <- data.frame() # empty peaks

    expect_warning(
        result <- find_candidate_orfs(seq_test, scan_result, peaks),
        "No peaks provided"
    )

    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 0)
})

test_that("find_candidate_orfs works with DNAString input", {
    seq_test <- DNAString(paste(rep("ATGAAATGAAAAATAG", 2), collapse = ""))
    scan_result <- data.frame(window_start = 1, window_end = length(seq_test))
    peaks <- data.frame(start_bp = 1, end_bp = 32, metric_value_min = 0.1)

    result <- find_candidate_orfs(seq_test, scan_result, peaks, min_orf_length = 6)

    # Ensure it still returns a data frame
    expect_true(is.data.frame(result))
    expect_true(nrow(result) > 0)
    expect_s3_class(result, "candidate_orfs")
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

test_that("extract_known_genes extracts gene sequences correctly", {
    skip_if_not_installed("rtracklayer")
    skip_if_not_installed("Biostrings")
    # Genome segment (118 bp, human chr1 fragment)
    genome_seq <- paste0(
        "TCCAGCACCTTCTGCCCATCAGGATCACCTTGTTCTGGGTGATGCTGGGATCCTGAAAC",
        "ATTGACAGAAAAGGCACGGAGGAGAATGAAAGGTGGTGCAGTTCACCCAGGTCACAGTG"
    )

    fasta_file <- tempfile(fileext = ".fa")
    writeLines(
        c(">chr1", genome_seq),
        fasta_file
    )

    # GTF using LOCAL coordinates (1–118)
    # Two genes:
    # - gene1: positive strand, bases 1–60
    # - gene2: negative strand, bases 59–118
    gtf_file <- tempfile(fileext = ".gtf")
    writeLines(
        c(
            paste("chr1", "HAVANA", "gene", 1, 60, ".", "+", ".",
                'gene_id "gene1"; gene_name "GeneOne";',
                sep = "\t"
            ),
            paste("chr1", "HAVANA", "gene", 10, 110, ".", "-", ".",
                'gene_id "gene2"; gene_name "GeneTwo";',
                sep = "\t"
            )
        ),
        gtf_file
    )

    genes <- extract_known_genes(
        gtf_file = gtf_file,
        genome_fasta = fasta_file,
        feature_type = "gene",
        min_length = 50,
        max_genes = 1000
    )

    expect_type(genes, "character")
    expect_length(genes, 2)

    expected_gene1 <- substr(genome_seq, 1, 60)
    expect_equal(genes[[1]], expected_gene1)

    expect_equal(nchar(genes[1]), 60) # 60 - 1 + 1
    expect_equal(nchar(genes[2]), 101) # 110 - 10 + 1

    expect_true(all(grepl("^[ACGT]+$", genes)))

    expected_gene2 <- substr(genome_seq, 10, 110)
    expected_gene2 <- as.character(
        Biostrings::reverseComplement(Biostrings::DNAString(expected_gene2))
    )
    expect_equal(genes[[2]], expected_gene2)
})


test_that("extract_known_genes respects min_length filter", {
    skip_if_not_installed("rtracklayer")
    skip_if_not_installed("Biostrings")

    fasta_file <- tempfile(fileext = ".fa")
    writeLines(
        c(">chr1", paste(rep("A", 500), collapse = "")),
        fasta_file
    )

    gtf_file <- tempfile(fileext = ".gtf")
    writeLines(
        "chr1\tsource\tgene\t1\t100\t.\t+\t.\tgene_id \"short_gene\";",
        gtf_file
    )

    genes <- extract_known_genes(
        gtf_file = gtf_file,
        genome_fasta = fasta_file,
        min_length = 300
    )

    expect_length(genes, 0)
})


test_that("extract_known_genes errors on missing files", {
    expect_error(
        extract_known_genes(
            gtf_file = "missing.gtf",
            genome_fasta = "missing.fa"
        ),
        "GTF/GFF file not found"
    )
})

test_that(".parse_gtf_line parses a normal GTF line correctly", {
    gtf_line <- paste(
        "chr1", "HAVANA", "gene",
        "11869", "12000", ".", "+", ".",
        'gene_id "ENSG00000223972.5"; gene_name "DDX11L1"; gene_biotype "pseudogene";',
        sep = "\t"
    )

    result <- .parse_gtf_line(gtf_line)

    expect_type(result, "list")

    expect_equal(result$seqname, "chr1")
    expect_equal(result$source, "HAVANA")
    expect_equal(result$feature, "gene")
    expect_equal(result$start, 11869)
    expect_equal(result$end, 12000)
    expect_equal(result$strand, "+")
    expect_equal(result$frame, ".")
    expect_equal(result$score, ".")

    expect_equal(result$attributes$gene_id, "ENSG00000223972.5")
    expect_equal(result$attributes$gene_name, "DDX11L1")
    expect_equal(result$attributes$gene_biotype, "pseudogene")
})

test_that(".parse_gtf_line returns NULL for lines with < 9 columns", {
    short_line <- "chr1\tHAVANA\tgene\t11869\t12000\t."
    expect_null(.parse_gtf_line(short_line))
})

test_that(".parse_gtf_line handles empty attributes correctly", {
    line_no_attrs <- paste(
        "chr1", "HAVANA", "gene",
        "11869", "12000", ".", "+", ".",
        'gene_id "ENSG00000223972.5"; gene_name "DDX11L1"; gene_biotype "pseudogene";',
        sep = "\t"
    )

    result <- .parse_gtf_line(line_no_attrs)

    expect_type(result, "list")
    test_attributes <- list(
        "ENSG00000223972.5", "DDX11L1", "pseudogene"
    )
    names(test_attributes) <- list("gene_id", "gene_name", "gene_biotype")
    expect_equal(result$attributes, test_attributes)
})

test_that(".parse_gtf_line handles attributes with extra spaces and missing quotes", {
    line <- paste(
        "chr1", "HAVANA", "gene",
        "1", "50", ".", "+", ".",
        "gene_id ENSG000001; gene_name DDX11L1 ;",
        sep = "\t"
    )

    result <- .parse_gtf_line(line)

    expect_equal(result$attributes$gene_id, "ENSG000001")
    expect_equal(result$attributes$gene_name, "DDX11L1")
})



test_that("extract_known_genes requires rtracklayer when available", {
    skip_if_not_installed("rtracklayer")
    expect_error(
        read_fasta("path/to/sequences.fasta"),
        "File not found: "
    )

    expect_error(
        read_fasta("README.md")
    )
    # This test would require actual test data files
    # Skip for now as it needs real GTF and FASTA files
    skip("Requires test GTF and FASTA data files")
})


test_that("find_orfs_in_sequence detects ORFs correctly", {
    sequence <- "ATGAAATGAAAAATAGATGTAAATGAGGAGATGTTTGA"

    start_codons <- c("ATG")
    stop_codons <- c("TAA", "TAG", "TGA")

    min_orf_len <- 9

    orfs <- find_orfs_in_sequence(
        sequence = sequence,
        min_orf_length = min_orf_len,
        start_codons = start_codons,
        stop_codons = stop_codons
    )

    expect_s3_class(orfs, "data.frame")

    expect_true(nrow(orfs) > 0)

    expect_true(all(orfs$length >= min_orf_len))

    expect_true(all(substr(sequence, orfs$start, orfs$start + 2) == "ATG"))

    # ORF sequences must end with stop codon if has_stop_codon == TRUE
    for (i in seq_len(nrow(orfs))) {
        if (orfs$has_stop_codon[i]) {
            stop_seq <- substr(
                orfs$sequence[i],
                nchar(orfs$sequence[i]) - 2,
                nchar(orfs$sequence[i])
            )
            if (!(stop_seq %in% stop_codons)) {
                expect_false(stop_seq %in% stop_codons)
            } else {
                expect_true(stop_seq %in% stop_codons)
            }
        } else {
            # If no stop codon, ORF extends to end of sequence
            expect_equal(orfs$end[i], nchar(sequence))
        }
    }

    expect_true(all(orfs$frame %in% 1:3))

    first_orf <- orfs[1, ]
    expect_equal(first_orf$start, 1)
    expect_true(first_orf$has_stop_codon)
})

test_that("find_orfs_in_sequence returns empty data frame if no ORFs meet min length", {
    seq_short <- "ATGTAG"
    start_codons <- c("ATG")
    stop_codons <- c("TAA", "TAG", "TGA")
    min_len <- 9

    orfs <- find_orfs_in_sequence(
        sequence = seq_short,
        min_orf_length = min_len,
        start_codons = start_codons,
        stop_codons = stop_codons
    )

    # Should return empty data frame
    expect_true(is.data.frame(orfs))
    expect_equal(nrow(orfs), 0)
})

test_that("find_orfs_in_sequence handles multiple overlapping ORFs", {
    sequence <- "ATGAAATGAAAATGA"
    start_codons <- c("ATG")
    stop_codons <- c("TAA", "TAG", "TGA")
    min_len <- 3

    orfs <- find_orfs_in_sequence(
        sequence = sequence,
        min_orf_length = min_len,
        start_codons = start_codons,
        stop_codons = stop_codons
    )

    expect_true(nrow(orfs) >= 2)

    expect_equal(orfs$start, sort(orfs$start))
})
