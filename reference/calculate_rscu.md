# Calculate Relative Synonymous Codon Usage (RSCU)

Calculate the Relative Synonymous Codon Usage for a codon frequency
distribution. RSCU \> 1 indicates codon usage bias towards that codon.

## Usage

``` r
calculate_rscu(frequencies)
```

## Arguments

- frequencies:

  Named numeric vector of codon frequencies

## Value

Named numeric vector of RSCU values

## Examples

``` r
sequence <- "ATGATGATGTTATTATTACGCCGCCGCC"
freqs <- calculate_codon_frequencies(sequence)
rscu <- calculate_rscu(freqs)
```
