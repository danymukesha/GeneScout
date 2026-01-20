# Compare Entropy Profiles

Compare entropy profiles between multiple sequences or regions.

## Usage

``` r
compare_entropy_profiles(..., metric = "shannon_entropy", labels = NULL)
```

## Arguments

- ...:

  One or more scan results or named lists with scan results

- metric:

  Character string: "shannon_entropy" or "kl_divergence"

- labels:

  Optional character vector of labels for each profile

## Value

A ggplot2 object

## Examples

``` r
seq1 <- paste(rep("ATGATGATG", 30), collapse = "")
seq2 <- paste(rep("ACGTACGT", 30), collapse = "")
result1 <- sliding_window_scan(seq1, window_size = 100, step_size = 20)
result2 <- sliding_window_scan(seq2, window_size = 100, step_size = 20)
plot <- compare_entropy_profiles(result1, result2, labels = c("Region 1", "Region 2"))
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the GeneScout package.
#>   Please report the issue at <https://github.com/danymukesha/GeneScout/issues>.
```
