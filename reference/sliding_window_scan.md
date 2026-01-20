# Sliding Window Entropy Scan

Perform a sliding window scan across a DNA sequence to calculate entropy
and KL divergence metrics at each position. This is the main function
for de novo ORF discovery.

## Usage

``` r
sliding_window_scan(
  sequence,
  window_size = 300,
  step_size = 30,
  reference_profile = NULL,
  min_codons = 10
)
```

## Arguments

- sequence:

  Character string or DNAString object representing DNA sequence

- window_size:

  Integer size of sliding window in base pairs (default: 300)

- step_size:

  Integer step size for sliding window (default: 30)

- reference_profile:

  Optional reference codon frequency profile. If NULL, only Shannon
  entropy is calculated.

- min_codons:

  Minimum number of complete codons required for analysis (default: 10)

## Value

A data frame with window coordinates and entropy metrics

## Examples

``` r
sequence <- paste(rep("ATGATGATGTTATTATTACGCCGCCGCC", 20), collapse = "")
result <- sliding_window_scan(sequence, window_size = 150, step_size = 30)
```
