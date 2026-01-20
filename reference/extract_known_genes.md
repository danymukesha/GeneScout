# Extract Known Genes from GTF/GFF File

Extract gene sequences from a GTF/GFF annotation file and corresponding
genome FASTA. This function parses the annotation file to get gene
coordinates and extracts sequences from the genome FASTA file.

## Usage

``` r
extract_known_genes(
  gtf_file,
  genome_fasta,
  feature_type = "gene",
  gene_attribute = "gene_id",
  min_length = 300,
  max_genes = Inf
)
```

## Arguments

- gtf_file:

  Character string specifying the path to the GTF/GFF file

- genome_fasta:

  Character string specifying the path to the genome FASTA file

- feature_type:

  Character string specifying the feature type to extract (default:
  "gene", alternatives: "exon", "CDS", "transcript")

- gene_attribute:

  Character string specifying the gene ID attribute (default: "gene_id")

- min_length:

  Minimum gene length in base pairs (default: 300)

- max_genes:

  Maximum number of genes to extract (default: Inf)

## Value

Character vector of gene sequences

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract gene sequences from GTF and genome files
known_genes <- extract_known_genes(
    gtf_file = "annotations.gtf",
    genome_fasta = "genome.fasta",
    min_length = 300
)

# Create reference profile
ref_profile <- create_reference_profile(known_genes)
} # }
```
