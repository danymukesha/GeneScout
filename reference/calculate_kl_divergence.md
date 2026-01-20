# Calculate Kullback-Leibler Divergence

Calculate the Kullback-Leibler (KL) divergence between two codon
frequency distributions. This measures how one distribution diverges
from a second, expected distribution.

## Usage

``` r
calculate_kl_divergence(observed, reference, epsilon = 1e-10)
```

## Arguments

- observed:

  Named numeric vector of observed codon frequencies

- reference:

  Named numeric vector of reference codon frequencies

- epsilon:

  Small value to avoid log(0) (default: 1e-10)

## Value

Numeric value of KL divergence

## Examples

``` r
sequence1 <- "ATGATGATGTTATTATTACGCCGCCGCC"
sequence2 <- "TTATTATTACGCCGCCGCCATGATGATG"
freqs1 <- calculate_codon_frequencies(sequence1)
freqs2 <- calculate_codon_frequencies(sequence2)
kl_div <- calculate_kl_divergence(freqs1, freqs2)
```
