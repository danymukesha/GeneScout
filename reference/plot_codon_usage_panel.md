# Plot Codon Usage Bias and RSCU in a Unified Panel

This function generates a clean, publication-ready multi-panel
visualization summarizing codon usage bias and relative synonymous codon
usage (RSCU). The figure integrates absolute usage, codon-level RSCU,
and amino-acid-level RSCU distributions while avoiding redundant
encodings.

## Usage

``` r
plot_codon_usage_panel(df)
```

## Arguments

- df:

  A data.frame or tibble with columns:

  - `codon`: character, codon sequence (e.g. "AAA")

  - `frequency`: numeric, codon usage frequency

  - `aa`: character, translated amino acid (single-letter code)

## Value

A patchwork object combining three ggplot panels.

## Details

Panels included:

1.  Absolute codon usage aggregated by amino acid

2.  Relative synonymous codon usage (RSCU) per codon, faceted by amino
    acid

3.  Distribution of RSCU values per amino acid

RSCU is computed as: \$\$RSCU =
\frac{f\_{codon}}{mean(f\_{synonymous})}\$\$

## Examples

``` r
if (FALSE) { # \dontrun{
plot_codon_usage_panel(df)
} # }
```
