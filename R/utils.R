#' Read FASTA File
#'
#' Read DNA sequences from a FASTA file and convert to Biostrings DNAStringSet.
#'
#' @param file_path Character string specifying the path to a FASTA or FASTA.GZ file
#' @return A \code{DNAStringSet} object containing the DNA sequences
#' @export
#' @examples
#' \dontrun{
#' sequences <- read_fasta("path/to/sequences.fasta")
#' }
read_fasta <- function(file_path) {
    if (!file.exists(file_path)) {
        stop("File not found: ", file_path)
    }

    dna_seqs <- Biostrings::readDNAStringSet(file_path)

    if (length(dna_seqs) == 0) {
        stop("No sequences found in FASTA file")
    }

    return(dna_seqs)
}

#' Get All Possible Codons
#'
#' Generate a vector of all 64 possible codons.
#'
#' @return Character vector of 64 codons
#' @keywords internal
.get_all_codons <- function() {
    bases <- c("A", "C", "G", "T")
    codons <- expand.grid(b1 = bases, b2 = bases, b3 = bases)
    return(paste(codons$b1, codons$b2, codons$b3, sep = ""))
}

#' Validate DNA Sequence
#'
#' Check if a DNA sequence contains only valid nucleotides.
#'
#' @param sequence Character string of DNA sequence
#' @return Boolean indicating if sequence is valid
#' @keywords internal
.validate_dna <- function(sequence) {
    valid_bases <- c("A", "C", "G", "T", "N")
    cleaned <- toupper(gsub("\\s", "", sequence))
    all_chars <- strsplit(cleaned, "")[[1]]
    return(all(all_chars %in% valid_bases))
}

#' Ensure Sequence Length is Divisible by 3
#'
#' Truncate or pad sequence to make length divisible by 3 for codon analysis.
#'
#' @param sequence Character string of DNA sequence
#' @param strategy Character string: "truncate" or "pad" (default: "truncate")
#' @return Character string of adjusted sequence
#' @keywords internal
.adjust_codon_length <- function(sequence, strategy = "truncate") {
    len <- nchar(sequence)
    remainder <- len %% 3

    if (remainder == 0) {
        return(sequence)
    }

    if (strategy == "truncate") {
        return(substr(sequence, 1, len - remainder))
    } else if (strategy == "pad") {
        padding <- paste(rep("N", 3 - remainder), collapse = "")
        return(paste0(sequence, padding))
    }

    stop("Invalid strategy. Use 'truncate' or 'pad'")
}

#' Sliding Window Iterator
#'
#' Helper function to generate sliding window coordinates.
#'
#' @param total_length Integer length of sequence
#' @param window_size Integer size of sliding window
#' @param step_size Integer step size for sliding window
#' @return Data frame with window coordinates
#' @keywords internal
.get_window_coords <- function(total_length, window_size, step_size) {
    num_windows <- floor((total_length - window_size) / step_size) + 1
    starts <- seq(1, total_length - window_size + 1, by = step_size)
    ends <- starts + window_size - 1
    return(data.frame(
        window_id = 1:length(starts),
        start = starts,
        end = ends,
        stringsAsFactors = FALSE
    ))
}

#' Extract Known Genes from GTF/GFF File
#'
#' Extract gene sequences from a GTF/GFF annotation file and corresponding genome FASTA.
#' This function parses the annotation file to get gene coordinates and extracts
#' sequences from the genome FASTA file.
#'
#' @param gtf_file Character string specifying the path to the GTF/GFF file
#' @param genome_fasta Character string specifying the path to the genome FASTA file
#' @param feature_type Character string specifying the feature type to extract
#'   (default: "gene", alternatives: "exon", "CDS", "transcript")
#' @param gene_attribute Character string specifying the gene ID attribute (default: "gene_id")
#' @param min_length Minimum gene length in base pairs (default: 300)
#' @param max_genes Maximum number of genes to extract (default: Inf)
#' @return Character vector of gene sequences
#' @export
#' @examples
#' \dontrun{
#' # Extract gene sequences from GTF and genome files
#' known_genes <- extract_known_genes(
#'     gtf_file = "annotations.gtf",
#'     genome_fasta = "genome.fasta",
#'     min_length = 300
#' )
#'
#' # Create reference profile
#' ref_profile <- create_reference_profile(known_genes)
#' }
extract_known_genes <- function(gtf_file,
                                genome_fasta,
                                feature_type = "gene",
                                gene_attribute = "gene_id",
                                min_length = 300,
                                max_genes = Inf) {
    if (!file.exists(gtf_file)) {
        stop("GTF/GFF file not found: ", gtf_file)
    }

    if (!file.exists(genome_fasta)) {
        stop("Genome FASTA file not found: ", genome_fasta)
    }

    if (!requireNamespace("rtracklayer", quietly = TRUE)) {
        stop(
            "Package 'rtracklayer' is required to parse GTF files. ",
            "Install with: BiocManager::install('rtracklayer')"
        )
    }

    gtf_data <- rtracklayer::import(gtf_file)

    genome <- Biostrings::readDNAStringSet(genome_fasta)
    genome_df <- data.frame(
        seqname = names(genome),
        seq = as.character(genome),
        stringsAsFactors = FALSE
    )

    if (feature_type %in% names(gtf_data)) {
        feature_col <- gtf_data[[feature_type]]
    } else {
        feature_col <- gtf_data$type
    }

    feature_rows <- which(feature_col == feature_type)

    if (length(feature_rows) == 0) {
        warning("No features of type '", feature_type, "' found in GTF file")
        return(character(0))
    }

    gene_sequences <- character(0)
    gene_count <- 0

    for (i in feature_rows) {
        if (gene_count >= max_genes) {
            break
        }

        chr <- as.character(gtf_data@seqnames[i])
        start <- as.numeric(gtf_data@ranges@start[i])
        end <- gtf_data@ranges@width[i] + gtf_data@ranges@start[i] - 1
        strand <- as.character(gtf_data@strand[i])

        genome_row <- which(sub(" .*", "", genome_df$seqname) == chr)

        if (length(genome_row) == 0) {
            next
        }

        genome_seq <- genome_df$seq[genome_row[1]]

        if (end > nchar(genome_seq)) {
            next
        }

        gene_seq <- substr(genome_seq, start, end)

        if (strand == "-") {
            gene_seq <- Biostrings::reverseComplement(Biostrings::DNAString(gene_seq))
            gene_seq <- as.character(gene_seq)
        }

        if (nchar(gene_seq) >= min_length) {
            if (.validate_dna(gene_seq)) {
                gene_sequences <- c(gene_sequences, gene_seq)
                gene_count <- gene_count + 1
            }
        }
    }

    if (length(gene_sequences) == 0) {
        warning("No gene sequences extracted. Check GTF file and parameters.")
    }

    message("Extracted ", length(gene_sequences), " gene sequences")

    return(gene_sequences)
}

#' Parse GTF Line
#'
#' Helper function to parse a single GTF line.
#'
#' @param line Character string of a GTF line
#' @return List of parsed attributes
#' @keywords internal
.parse_gtf_line <- function(line) {
    fields <- strsplit(line, "\t")[[1]]

    if (length(fields) < 9) {
        return(NULL)
    }

    seqname <- fields[1]
    source <- fields[2]
    feature <- fields[3]
    start <- as.numeric(fields[4])
    end <- as.numeric(fields[5])
    score <- fields[6]
    strand <- fields[7]
    frame <- fields[8]
    attributes <- fields[9]

    attr_dict <- list()
    attr_pairs <- strsplit(attributes, ";")[[1]]
    attr_pairs <- trimws(attr_pairs)
    attr_pairs <- attr_pairs[attr_pairs != ""]

    for (pair in attr_pairs) {
        key_value <- strsplit(pair, " ", fixed = TRUE)[[1]]
        if (length(key_value) >= 2) {
            key <- key_value[1]
            value <- paste(key_value[-1], collapse = " ")
            value <- gsub('"', "", value)
            attr_dict[[key]] <- value
        }
    }

    return(list(
        seqname = seqname,
        source = source,
        feature = feature,
        start = start,
        end = end,
        score = score,
        strand = strand,
        frame = frame,
        attributes = attr_dict
    ))
}
