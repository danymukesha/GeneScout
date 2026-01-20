# Calculate Effective Number of Codons (ENC)

Calculate the Effective Number of Codons (ENC), a measure of codon usage
bias. ENC ranges from 20 (extreme bias) to 61 (no bias).

## Usage

``` r
calculate_enc(frequencies)
```

## Arguments

- frequencies:

  Named numeric vector of codon frequencies

## Value

Numeric value of ENC

## Examples

``` r
sequence <- "ATGATGATGTTATTATTACGCCGCCGCC"
freqs <- calculate_codon_frequencies(sequence)
enc <- calculate_enc(freqs)
#> [1] 6
#> [1] "enc: 0 entro: 0"
#> [1] 6
#> [1] "enc: 0 entro: 0"
#> [1] 1
#> [1] "enc: 0 entro: 0"
```
