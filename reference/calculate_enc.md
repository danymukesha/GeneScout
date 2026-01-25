# Calculate the Effective Number of Codons (ENC)

Computes the Effective Number of Codons (ENC; Wright 1990), a widely
used measure of synonymous codon usage bias in protein-coding DNA
sequences. ENC quantifies how evenly synonymous codons are used across
amino acids.

## Usage

``` r
calculate_enc(frequencies)
```

## Arguments

- frequencies:

  Named numeric vector of codon frequencies or counts. Names must be
  standard DNA codons (e.g., "ATG", "GCC") and should include all 64
  codons. Values may be raw counts or normalized frequencies.

## Value

A single numeric value giving the Effective Number of Codons (ENC),
constrained to a maximum of 61.

## Details

ENC ranges from:

- **20**: extreme codon usage bias (one codon per amino acid)

- **61**: no codon usage bias (all synonymous codons used equally)

This implementation follows the original formulation by Wright (1990),
based on codon-family homozygosity (\\F_k\\) for amino acids with \\k =
2, 3, 4, 6\\ synonymous codons:

\$\$ ENC = 2 + \frac{9}{F_2} + \frac{1}{F_3} + \frac{5}{F_4} +
\frac{3}{F_6} \$\$

where \\F_k\\ is the average homozygosity of codon usage within each
synonymous codon family of size \\k\\.

Stop codons and amino acids encoded by a single codon (Methionine and
Tryptophan) are excluded from the calculation, as they do not contribute
to synonymous codon bias.

The input vector is internally normalized within each synonymous codon
family. Amino acids for which no codons are observed are ignored. If
insufficient information is available for one or more degeneracy
classes, their contribution to ENC is omitted.

## References

Wright, F. (1990). The 'effective number of codons' used in a gene.
*Gene*, 87(1), 23–29. doi:10.1016/0378-1119(90)90491-9

## See also

[`calculate_codon_frequencies`](https://danymukesha.github.io/GeneScout/reference/calculate_codon_frequencies.md),
[`calculate_rscu`](https://danymukesha.github.io/GeneScout/reference/calculate_rscu.md)

## Examples

``` r
sequence <- "ATGATGATGTTATTATTACGCCGCCGCC"
freqs <- calculate_codon_frequencies(sequence)
calculate_enc(freqs)
#> Warning: argument is not numeric or logical: returning NA
#> Warning: argument is not numeric or logical: returning NA
#> Warning: argument is not numeric or logical: returning NA
#> [1] 5
```
