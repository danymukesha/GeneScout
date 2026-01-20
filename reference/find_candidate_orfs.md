# Find Candidate ORFs

Identify candidate open reading frames (ORFs) based on entropy analysis
and sequence features.

## Usage

``` r
find_candidate_orfs(
  sequence,
  scan_result,
  peaks,
  min_orf_length = 90,
  start_codons = c("ATG"),
  stop_codons = c("TAA", "TAG", "TGA")
)
```

## Arguments

- sequence:

  Character string or DNAString object representing DNA sequence

- scan_result:

  Data frame from sliding_window_scan()

- peaks:

  Data frame from entropy_peak_detection()

- min_orf_length:

  Minimum ORF length in base pairs (default: 90, i.e., 30 aa)

- start_codons:

  Character vector of valid start codons (default: c("ATG"))

- stop_codons:

  Character vector of stop codons (default: c("TAA", "TAG", "TGA"))

## Value

Data frame with candidate ORF information

## Examples

``` r
sequence <- paste(rep("ATGATGATGTTATTATTACGCCGCCGCC", 20), collapse = "")
result <- sliding_window_scan(sequence, window_size = 150, step_size = 30)
peaks <- entropy_peak_detection(result, threshold = 0.1)
#> Warning: No peaks detected with current threshold
candidates <- find_candidate_orfs(sequence, result, peaks)
#> Warning: No peaks provided, returning empty data frame
```
