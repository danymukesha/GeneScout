# Ensure Sequence Length is Divisible by 3

Truncate or pad sequence to make length divisible by 3 for codon
analysis.

## Usage

``` r
.adjust_codon_length(sequence, strategy = "truncate")
```

## Arguments

- sequence:

  Character string of DNA sequence

- strategy:

  Character string: "truncate" or "pad" (default: "truncate")

## Value

Character string of adjusted sequence
