# Calculate Codon Frequencies

Calculate the frequency distribution of all 64 codons in a DNA sequence.

## Usage

``` r
calculate_codon_frequencies(sequence, normalize = TRUE)
```

## Arguments

- sequence:

  Character string or DNAString object representing DNA sequence

- normalize:

  Boolean indicating whether to normalize frequencies (default: TRUE)

## Value

Named numeric vector of codon frequencies

## Examples

``` r
sequence <- "ATGATGATGTTATTATTACGCCGCCGCC"
frequencies <- calculate_codon_frequencies(sequence)
```
