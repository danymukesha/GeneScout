# Calculate Shannon Entropy

Calculate Shannon entropy of a codon frequency distribution.

## Usage

``` r
calculate_shannon_entropy(frequencies)
```

## Arguments

- frequencies:

  Named numeric vector of codon frequencies

## Value

Numeric value of Shannon entropy in bits

## Examples

``` r
sequence <- "ATGATGATGTTATTATTACGCCGCCGCC"
freqs <- calculate_codon_frequencies(sequence)
entropy <- calculate_shannon_entropy(freqs)
```
