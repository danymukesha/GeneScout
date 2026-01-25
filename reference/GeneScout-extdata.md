# Example FASTA and GTF files bundled with GeneScout

GeneScout includes a collection of small to moderately sized FASTA and
GTF files intended for examples, testing, and vignettes. These files are
stored in compressed form under `inst/extdata` to reduce package size
and improve portability.

## Details

The datasets include:

- Human protein-coding and lncRNA transcript FASTA files from GENCODE
  v19

- Human chromosome 19 genomic sequence (GRCh37) \[`REMOVED CURRENTLY!`\]

- Example gene sequences (e.g., `APOE`, `TP53`)

- A nematode genome and gene annotation (Acanthocheilonema viteae)

Files can be accessed using
[`system.file`](https://rdrr.io/r/base/system.file.html):

    system.file("extdata", "APOE.fasta.gz", package = "GeneScout")

All files are provided strictly for demonstration and benchmarking
purposes and should not be considered authoritative reference genomes.
