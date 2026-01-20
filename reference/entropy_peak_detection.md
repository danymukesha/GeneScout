# Detect Entropy Peaks

Identify peaks in entropy profiles that may indicate potential coding
regions.

## Usage

``` r
entropy_peak_detection(
  scan_result,
  metric = "shannon_entropy",
  method = "quantile",
  threshold = 0.1,
  min_peak_width = 3
)
```

## Arguments

- scan_result:

  Data frame from sliding_window_scan()

- metric:

  Character string: "shannon_entropy" or "kl_divergence"

- method:

  Character string: "quantile" or "sd" for threshold calculation

- threshold:

  Quantile threshold (for "quantile" method) or number of SDs (for "sd"
  method)

- min_peak_width:

  Minimum width of peaks in windows (default: 3)

## Value

Data frame with peak information

## Examples

``` r
sequence <- paste(rep("ATGATGATGTTATTATTACGCCGCCGCC", 20), collapse = "")
result <- sliding_window_scan(sequence, window_size = 150, step_size = 30)
peaks <- entropy_peak_detection(result, threshold = 0.1)
#> Warning: No peaks detected with current threshold
```
