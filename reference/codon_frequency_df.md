# Convert Codon Frequency Vector to Tidy Data Frame

Convert Codon Frequency Vector to Tidy Data Frame

## Usage

``` r
codon_frequency_df(frequencies)
```

## Arguments

- frequencies:

  Named numeric vector of codon frequencies

## Value

Tibble with columns codon, frequency, aa

## Examples

``` r
freqs <- calculate_codon_frequencies("ATGATGATG")
codon_frequency_df(freqs)
#>    codon frequency aa
#> 1    AAA         0  K
#> 2    CAA         0  Q
#> 3    GAA         0  E
#> 4    TAA         0  *
#> 5    ACA         0  T
#> 6    CCA         0  P
#> 7    GCA         0  A
#> 8    TCA         0  S
#> 9    AGA         0  R
#> 10   CGA         0  R
#> 11   GGA         0  G
#> 12   TGA         0  *
#> 13   ATA         0  I
#> 14   CTA         0  L
#> 15   GTA         0  V
#> 16   TTA         0  L
#> 17   AAC         0  N
#> 18   CAC         0  H
#> 19   GAC         0  D
#> 20   TAC         0  Y
#> 21   ACC         0  T
#> 22   CCC         0  P
#> 23   GCC         0  A
#> 24   TCC         0  S
#> 25   AGC         0  S
#> 26   CGC         0  R
#> 27   GGC         0  G
#> 28   TGC         0  C
#> 29   ATC         0  I
#> 30   CTC         0  L
#> 31   GTC         0  V
#> 32   TTC         0  F
#> 33   AAG         0  K
#> 34   CAG         0  Q
#> 35   GAG         0  E
#> 36   TAG         0  *
#> 37   ACG         0  T
#> 38   CCG         0  P
#> 39   GCG         0  A
#> 40   TCG         0  S
#> 41   AGG         0  R
#> 42   CGG         0  R
#> 43   GGG         0  G
#> 44   TGG         0  W
#> 45   ATG         1  M
#> 46   CTG         0  L
#> 47   GTG         0  V
#> 48   TTG         0  L
#> 49   AAT         0  N
#> 50   CAT         0  H
#> 51   GAT         0  D
#> 52   TAT         0  Y
#> 53   ACT         0  T
#> 54   CCT         0  P
#> 55   GCT         0  A
#> 56   TCT         0  S
#> 57   AGT         0  S
#> 58   CGT         0  R
#> 59   GGT         0  G
#> 60   TGT         0  C
#> 61   ATT         0  I
#> 62   CTT         0  L
#> 63   GTT         0  V
#> 64   TTT         0  F
```
