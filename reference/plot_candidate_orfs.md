# Plot Candidate ORFs

Visualize candidate ORFs along with entropy profile.

## Usage

``` r
plot_candidate_orfs(
  scan_result,
  candidates,
  peaks = NULL,
  metric = "shannon_entropy"
)
```

## Arguments

- scan_result:

  Data frame from sliding_window_scan()

- candidates:

  Data frame from find_candidate_orfs()

- peaks:

  Optional data frame from entropy_peak_detection()

- metric:

  Character string: "shannon_entropy" or "kl_divergence"

## Value

A ggplot2 object

## Examples

``` r
sequence <- paste(rep("ATGATGATGTTATTATTACGCCGCCGCC", 20), collapse = "")
result <- sliding_window_scan(sequence, window_size = 150, step_size = 30)
peaks <- entropy_peak_detection(result, threshold = 0.1)
#> Warning: No peaks detected with current threshold
candidates <- find_candidate_orfs(sequence, result, peaks)
#> Warning: No peaks provided, returning empty data frame
plot <- plot_candidate_orfs(result, candidates, peaks)
```
