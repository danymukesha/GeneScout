#' Example FASTA and GTF files bundled with GeneScout
#'
#' GeneScout includes a collection of small to moderately sized FASTA and GTF
#' files intended for examples, testing, and vignettes. These files are stored
#' in compressed form under \code{inst/extdata} to reduce package size and
#' improve portability.
#'
#' The datasets include:
#' \itemize{
#'   \item Human protein-coding and lncRNA transcript FASTA files from GENCODE v19
#'   \item Human chromosome 19 genomic sequence (GRCh37) \[`REMOVED CURRENTLY!`\]
#'   \item Example gene sequences (e.g., \code{APOE}, \code{TP53})
#'   \item A nematode genome and gene annotation (Acanthocheilonema viteae)
#' }
#'
#' Files can be accessed using \code{\link[base]{system.file}}:
#'
#' \preformatted{
#' system.file("extdata", "APOE.fasta.gz", package = "GeneScout")
#' }
#'
#' All files are provided strictly for demonstration and benchmarking purposes
#' and should not be considered authoritative reference genomes.
#'
#' @name GeneScout-extdata
#' @docType data
#' @keywords datasets
NULL
