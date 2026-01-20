# Create Reference Codon Profile

Create a reference codon usage profile from known genes or a set of
training sequences. This profile represents the organism's codon usage
bias pattern.

## Usage

``` r
create_reference_profile(sequences, method = "mean")
```

## Arguments

- sequences:

  Character vector, DNAStringSet, or list of DNA sequences

- method:

  Character string: "mean" or "median" for aggregating frequencies

## Value

Named numeric vector representing reference codon frequencies

## Examples

``` r
known_genes <- c("ATGATGATG", "GCCGCCGCC", "TTATTATTA")
ref_profile <- create_reference_profile(known_genes)
```
