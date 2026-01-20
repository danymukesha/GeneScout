
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GeneScout

<!-- badges: start -->

<!-- badges: end -->

**GeneScout** is an R package for identifying hidden genes in non-coding
DNA using statistical information theory. It implements sliding window
entropy scanning, Shannon entropy calculation, and Kullback-Leibler
divergence to discover potential small open reading frames (ORFs) in
genomic sequences without requiring prior annotation.

Most existing bioinformatics tools for codon usage bias analysis (like
`coRdon` and `Biostrings`) are designed for analyzing **known** genes
from annotated GTF/GFF files. They are not optimized for **discovery**
in non-coding regions.

**GeneScout fills this gap** by focusing on *de novo* small ORF
discovery in large genomic regions.

## Features

- **Sliding Window Entropy Scanning**: Efficiently scan millions of base
  pairs to find regions with unusual codon usage patterns
- **Shannon Entropy Calculation**: Measure the randomness of codon
  distributions
- **Kullback-Leibler Divergence**: Compare regions to organism-specific
  codon usage profiles
- **ORF Detection**: Find candidate open reading frames in low-entropy
  regions
- **Reference Profile Creation**: Build organism-specific codon usage
  profiles from known genes

## Installation

``` r
# install development version
devtools::install_github("danymukesha/GeneScout")
```

## Quick Start

``` r
library(GeneScout)

# Create a reference profile from known genes
known_genes <- c("ATGATGATG", "GCCGCCGCC", "TTATTATTA")
ref_profile <- create_reference_profile(known_genes)

# Scan a genomic sequence
sequence <- paste(rep("ATGATGATGTTATTATTACGCCGCCGCC", 30), collapse = "")
scan_result <- sliding_window_scan(
  sequence,
  window_size = 300,
  step_size = 30,
  reference_profile = ref_profile
)

# Detect entropy peaks
peaks <- entropy_peak_detection(scan_result, threshold = 0.1)

# Find candidate ORFs
candidates <- find_candidate_orfs(sequence, scan_result, peaks)

# Visualize results
plot_candidate_orfs(scan_result, candidates, peaks)
```

<!-- ## Core Functions -->

<!-- ### Entropy Analysis -->

<!-- - `calculate_codon_frequencies()` - Calculate codon usage frequencies -->

<!-- - `calculate_shannon_entropy()` - Calculate Shannon entropy -->

<!-- - `calculate_kl_divergence()` - Calculate KL divergence -->

<!-- - `calculate_rscu()` - Calculate Relative Synonymous Codon Usage -->

<!-- - `calculate_enc()` - Calculate Effective Number of Codons -->

<!-- ### Discovery -->

<!-- - `read_fasta()` - Read DNA sequences from FASTA files -->

<!-- - `extract_known_genes()` - Extract gene sequences from GTF/GFF annotations -->

<!-- - `create_reference_profile()` - Create organism-specific reference profile -->

<!-- - `sliding_window_scan()` - Main scanning function -->

<!-- - `entropy_peak_detection()` - Detect entropy peaks -->

<!-- - `find_candidate_orfs()` - Find candidate ORFs -->

<!-- ### Visualization -->

<!-- - `plot_entropy_profile()` - Plot entropy profiles -->

<!-- - `compare_entropy_profiles()` - Compare multiple profiles -->

<!-- - `plot_candidate_orfs()` - Visualize candidate ORFs -->

## Theory

### Shannon Entropy

$$H = -\sum_{i=1}^{64} P(i) \log_2 P(i)$$

- **High entropy** (~6 bits): Random codon usage → non-coding DNA
- **Low entropy** (~3-4 bits): Biased codon usage → potential gene

### Kullback-Leibler Divergence

$$D_{KL}(P \| Q) = \sum_{i} P(i) \log_2 \left( \frac{P(i)}{Q(i)} \right)$$

Compares the observed codon distribution ($P$) to the reference
distribution ($Q$).

## Use Cases

1.  **De Novo Gene Discovery**: Find potential coding regions in
    unannotated genomes
2.  **Small ORF Identification**: Discover small peptides in non-coding
    regions
3.  **Comparative Genomics**: Compare codon usage between organisms
4.  **RNA-Seq Integration**: Validate candidate ORFs with expression
    data
5.  **Evolutionary Studies**: Analyze codon usage patterns across
    lineages

## Acknowledgments

Built on top of Bioconductor’s excellent `Biostrings`, `coRdon`, and
`ggplot2` packages.
