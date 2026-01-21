# Changelog

## Version 0.1.0 (2026-01-20)

NEW FEATURES

- Added calculate_codon_frequencies() for codon usage analysis
- Added calculate_shannon_entropy() for entropy calculation
- Added calculate_kl_divergence() for comparing codon distributions
- Added calculate_rscu() for Relative Synonymous Codon Usage
- Added calculate_enc() for Effective Number of Codons
- Added read_fasta() for reading FASTA files
- Added extract_known_genes() for extracting genes from GTF/GFF
  annotations
- Added create_reference_profile() for building organism-specific
  profiles
- Added sliding_window_scan() for genome-wide entropy scanning
- Added entropy_peak_detection() for identifying low-entropy regions
- Added find_candidate_orfs() for ORF discovery
- Added plot_entropy_profile() for visualization
- Added compare_entropy_profiles() for comparative analysis
- Added plot_candidate_orfs() for candidate ORF visualization

DOCUMENTATION

- Added comprehensive vignette with examples
- Added function documentation with roxygen2
- Added example datasets (chr19.fasta, APOE.fasta)

TESTS

- Added testthat test suite with 20+ tests
- Tests cover all major functions

BUG FIXES

- None (initial release)

DEPRECATED AND DEFUNCT

- None
