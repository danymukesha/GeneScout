# Find ORFs in Sequence

Helper function to find ORFs within a specific sequence region.

## Usage

``` r
find_orfs_in_sequence(sequence, min_orf_length, start_codons, stop_codons)
```

## Arguments

- sequence:

  Character string of DNA sequence

- min_orf_length:

  Minimum ORF length in base pairs

- start_codons:

  Character vector of start codons

- stop_codons:

  Character vector of stop codons

## Value

Data frame with ORF information
