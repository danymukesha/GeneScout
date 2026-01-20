#' Plot Entropy Profile
#'
#' Create a visualization of entropy metrics across the scanned sequence.
#'
#' @param scan_result Data frame from sliding_window_scan()
#' @param peaks Optional data frame from entropy_peak_detection()
#' @param metric Character string: "shannon_entropy" or "kl_divergence"
#' @param highlight_threshold Boolean indicating whether to show threshold line
#' @param show_peaks Boolean indicating whether to highlight detected peaks
#' @return A ggplot2 object
#' @export
#' @examples
#' sequence <- paste(rep("ATGATGATGTTATTATTACGCCGCCGCC", 20), collapse = "")
#' result <- sliding_window_scan(sequence, window_size = 150, step_size = 30)
#' plot <- plot_entropy_profile(result)
plot_entropy_profile <- function(scan_result,
                                   peaks = NULL,
                                   metric = "shannon_entropy",
                                   highlight_threshold = TRUE,
                                   show_peaks = FALSE) {
  if (!inherits(scan_result, "entropy_scan_result")) {
    scan_result <- structure(scan_result, class = c("entropy_scan_result", class(scan_result)))
  }

  if (!metric %in% c("shannon_entropy", "kl_divergence")) {
    stop("Invalid metric. Use 'shannon_entropy' or 'kl_divergence'")
  }

  if (metric == "kl_divergence" && all(is.na(scan_result$kl_divergence))) {
    stop("KL divergence not available. Provide reference_profile to sliding_window_scan()")
  }

  plot_data <- dplyr::mutate(scan_result,
    mid_position = (start + end) / 2
  )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = mid_position, y = .data[[metric]])) +
    ggplot2::geom_line(alpha = 0.7, color = "#2E86AB") +
    ggplot2::geom_point(alpha = 0.4, size = 0.8, color = "#2E86AB") +
    ggplot2::labs(
      x = "Position (bp)",
      y = if (metric == "shannon_entropy") "Shannon Entropy (bits)" else "KL Divergence",
      title = paste("Entropy Profile -", toupper(gsub("_", " ", metric)))
    ) +
    ggplot2::theme_minimal()

  if (highlight_threshold) {
    if (metric == "shannon_entropy") {
      cutoff <- quantile(scan_result[[metric]], 0.1, na.rm = TRUE)
      p <- p + ggplot2::geom_hline(yintercept = cutoff, 
                                     linetype = "dashed", 
                                     color = "#A23B72",
                                     alpha = 0.6)
    } else {
      cutoff <- quantile(scan_result[[metric]], 0.9, na.rm = TRUE)
      p <- p + ggplot2::geom_hline(yintercept = cutoff, 
                                     linetype = "dashed", 
                                     color = "#A23B72",
                                     alpha = 0.6)
    }
  }

  if (show_peaks && !is.null(peaks) && nrow(peaks) > 0) {
    for (i in 1:nrow(peaks)) {
      p <- p + ggplot2::geom_rect(
        ggplot2::aes(xmin = peaks$start_bp[i], 
                      xmax = peaks$end_bp[i],
                      ymin = -Inf, 
                      ymax = Inf),
        fill = "#F18F01",
        alpha = 0.2,
        inherit.aes = FALSE
      )
    }
  }

  return(p)
}

#' @export
plot.entropy_scan_result <- function(x, y, ...) {
  plot_entropy_profile(x, ...)
}

#' Compare Entropy Profiles
#'
#' Compare entropy profiles between multiple sequences or regions.
#'
#' @param ... One or more scan results or named lists with scan results
#' @param metric Character string: "shannon_entropy" or "kl_divergence"
#' @param labels Optional character vector of labels for each profile
#' @return A ggplot2 object
#' @export
#' @examples
#' seq1 <- paste(rep("ATGATGATG", 30), collapse = "")
#' seq2 <- paste(rep("ACGTACGT", 30), collapse = "")
#' result1 <- sliding_window_scan(seq1, window_size = 100, step_size = 20)
#' result2 <- sliding_window_scan(seq2, window_size = 100, step_size = 20)
#' plot <- compare_entropy_profiles(result1, result2, labels = c("Region 1", "Region 2"))
compare_entropy_profiles <- function(..., metric = "shannon_entropy", labels = NULL) {
  scan_results <- list(...)

  if (is.null(labels)) {
    labels <- paste("Profile", 1:length(scan_results))
  }

  if (length(labels) != length(scan_results)) {
    stop("Number of labels must match number of scan results")
  }

  plot_data <- lapply(1:length(scan_results), function(i) {
    df <- scan_results[[i]]
    df <- dplyr::mutate(df,
      mid_position = (start + end) / 2,
      profile = labels[i]
    )
    return(df)
  })

  combined_data <- dplyr::bind_rows(plot_data)

  colors <- c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D", "#3B1F2B")[1:length(labels)]

  p <- ggplot2::ggplot(combined_data, 
                        ggplot2::aes(x = mid_position, 
                                      y = .data[[metric]], 
                                      color = profile)) +
    ggplot2::geom_line(alpha = 0.7, size = 0.8) +
    ggplot2::geom_point(alpha = 0.3, size = 0.6) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(
      x = "Position (bp)",
      y = if (metric == "shannon_entropy") "Shannon Entropy (bits)" else "KL Divergence",
      title = paste("Comparison of", toupper(gsub("_", " ", metric)))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  return(p)
}

#' Plot Candidate ORFs
#'
#' Visualize candidate ORFs along with entropy profile.
#'
#' @param scan_result Data frame from sliding_window_scan()
#' @param candidates Data frame from find_candidate_orfs()
#' @param peaks Optional data frame from entropy_peak_detection()
#' @param metric Character string: "shannon_entropy" or "kl_divergence"
#' @return A ggplot2 object
#' @export
#' @examples
#' sequence <- paste(rep("ATGATGATGTTATTATTACGCCGCCGCC", 20), collapse = "")
#' result <- sliding_window_scan(sequence, window_size = 150, step_size = 30)
#' peaks <- entropy_peak_detection(result, threshold = 0.1)
#' candidates <- find_candidate_orfs(sequence, result, peaks)
#' plot <- plot_candidate_orfs(result, candidates, peaks)
plot_candidate_orfs <- function(scan_result,
                                 candidates,
                                 peaks = NULL,
                                 metric = "shannon_entropy") {
  if (!inherits(scan_result, "entropy_scan_result")) {
    scan_result <- structure(scan_result, class = c("entropy_scan_result", class(scan_result)))
  }

  plot_data <- dplyr::mutate(scan_result,
    mid_position = (start + end) / 2
  )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = mid_position, y = .data[[metric]])) +
    ggplot2::geom_line(alpha = 0.7, color = "#2E86AB") +
    ggplot2::labs(
      x = "Position (bp)",
      y = if (metric == "shannon_entropy") "Shannon Entropy (bits)" else "KL Divergence",
      title = "Candidate ORFs"
    ) +
    ggplot2::theme_minimal()

  if (!is.null(peaks) && nrow(peaks) > 0) {
    for (i in 1:nrow(peaks)) {
      p <- p + ggplot2::geom_rect(
        ggplot2::aes(xmin = peaks$start_bp[i], 
                      xmax = peaks$end_bp[i],
                      ymin = -Inf, 
                      ymax = Inf),
        fill = "#F18F01",
        alpha = 0.15,
        inherit.aes = FALSE
      )
    }
  }

  if (!is.null(candidates) && nrow(candidates) > 0) {
    orf_data <- dplyr::mutate(candidates,
      orf_id = 1:nrow(candidates)
    )

    p <- p + ggplot2::geom_rect(
      data = orf_data,
      ggplot2::aes(xmin = start, 
                    xmax = end,
                    ymin = -Inf, 
                    ymax = Inf,
                    group = orf_id),
      fill = "#A23B72",
      alpha = 0.3,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = orf_data,
      ggplot2::aes(x = (start + end) / 2,
                    y = Inf,
                    label = paste0("ORF", orf_id)),
      vjust = 1.5,
      hjust = 0.5,
      size = 3,
      color = "#A23B72",
      inherit.aes = FALSE
    )
  }

  return(p)
}
