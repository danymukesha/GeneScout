# Plot Entropy Profile

Create a visualization of entropy metrics across the scanned sequence.

## Usage

``` r
plot_entropy_profile(
  scan_result,
  peaks = NULL,
  metric = "shannon_entropy",
  highlight_threshold = TRUE,
  show_peaks = FALSE
)
```

## Arguments

- scan_result:

  Data frame from sliding_window_scan()

- peaks:

  Optional data frame from entropy_peak_detection()

- metric:

  Character string: "shannon_entropy" or "kl_divergence"

- highlight_threshold:

  Boolean indicating whether to show threshold line

- show_peaks:

  Boolean indicating whether to highlight detected peaks

## Value

A ggplot2 object

## Examples

``` r
sequence <- paste(rep("ATGATGATGTTATTATTACGCCGCCGCC", 20), collapse = "")
result <- sliding_window_scan(sequence, window_size = 150, step_size = 30)
plot <- plot_entropy_profile(result)
```
