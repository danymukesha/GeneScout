# Translate Codons to Amino Acids

Translate DNA codons into single-letter amino acid codes using the
standard genetic code.

## Usage

``` r
translate_codons(codons)
```

## Arguments

- codons:

  Character vector of codons (e.g. "ATG", "GCC")

## Value

Character vector of amino acids (single-letter codes, "\*" for stop)

## Examples

``` r
translate_codons(c("ATG", "TAA", "GCC"))
#> [1] "M" "*" "A"
```
