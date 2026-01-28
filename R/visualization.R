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
            p <- p + ggplot2::geom_hline(
                yintercept = cutoff,
                linetype = "dashed",
                color = "#A23B72",
                alpha = 0.6
            )
        } else {
            cutoff <- quantile(scan_result[[metric]], 0.9, na.rm = TRUE)
            p <- p + ggplot2::geom_hline(
                yintercept = cutoff,
                linetype = "dashed",
                color = "#A23B72",
                alpha = 0.6
            )
        }
    }

    if (show_peaks && !is.null(peaks) && nrow(peaks) > 0) {
        for (i in 1:nrow(peaks)) {
            p <- p + ggplot2::geom_rect(
                ggplot2::aes(
                    xmin = peaks$start_bp[i],
                    xmax = peaks$end_bp[i],
                    ymin = -Inf,
                    ymax = Inf
                ),
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

    p <- ggplot2::ggplot(
        combined_data,
        ggplot2::aes(
            x = mid_position,
            y = .data[[metric]],
            color = profile
        )
    ) +
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
            p <- p + ggplot2::annotate(
                "rect",
                xmin = peaks$start_bp[i],
                xmax = peaks$end_bp[i],
                ymin = -Inf,
                ymax = Inf,
                alpha = 0.15,
                fill = "#F18F01"
            )
        }
    }

    if (!is.null(candidates) && nrow(candidates) > 0) {
        orf_data <- dplyr::mutate(candidates,
            orf_id = 1:nrow(candidates),
            mid_position = (start + end) / 2
        )

        p <- p +
            ggplot2::geom_rect(
                data = orf_data,
                ggplot2::aes(
                    xmin = start,
                    xmax = end,
                    ymin = -Inf,
                    ymax = Inf,
                    group = orf_id
                ),
                fill = "#A23B72",
                alpha = 0.3,
                inherit.aes = FALSE
            ) +
            ggrepel::geom_text_repel(
                data = orf_data,
                ggplot2::aes(
                    x = mid_position,
                    y = max(plot_data[[metric]], na.rm = TRUE), # place above the plot
                    label = paste0("ORF", orf_id)
                ),
                nudge_y = 0.05 * max(plot_data[[metric]], na.rm = TRUE), # small offset
                direction = "x",
                angle = 0,
                vjust = 0,
                hjust = 0.5,
                size = 3,
                color = "black",
                inherit.aes = FALSE,
                segment.color = "grey30" # optional line connecting label to ORF
            )
    }

    return(p)
}

#' Plot Codon Usage Bias and RSCU in a Unified Panel
#'
#' This function generates a clean, publication-ready multi-panel visualization
#' summarizing codon usage bias and relative synonymous codon usage (RSCU).
#' The figure integrates absolute usage, codon-level RSCU, and amino-acid-level
#' RSCU distributions while avoiding redundant encodings.
#'
#' @param df A data.frame or tibble with columns:
#'   \itemize{
#'     \item \code{codon}: character, codon sequence (e.g. "AAA")
#'     \item \code{frequency}: numeric, codon usage frequency
#'     \item \code{aa}: character, translated amino acid (single-letter code)
#'   }
#'
#' @return A patchwork object combining three ggplot panels.
#'
#' @details
#' Panels included:
#' \enumerate{
#'   \item Absolute codon usage aggregated by amino acid
#'   \item Relative synonymous codon usage (RSCU) per codon, faceted by amino acid
#'   \item Distribution of RSCU values per amino acid
#' }
#'
#' RSCU is computed as:
#' \deqn{RSCU = \frac{f_{codon}}{mean(f_{synonymous})}}
#'
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#'
#' @examples
#' \dontrun{
#' plot_codon_usage_panel(df)
#' }
#'
#' @export
plot_codon_usage_panel <- function(df) {
    stopifnot(
        is.data.frame(df),
        all(c("codon", "frequency", "aa") %in% colnames(df)),
        is.numeric(df$frequency),
        all(df$frequency >= 0)
    )

    df <- df |>
        dplyr::mutate(
            codon = factor(codon),
            aa    = factor(aa)
        )

    df_rscu <- df |>
        dplyr::group_by(aa) |>
        dplyr::mutate(
            rscu = frequency / mean(frequency)
        ) |>
        dplyr::ungroup()

    p1 <- ggplot2::ggplot(df, ggplot2::aes(x = aa, y = frequency)) +
        ggplot2::geom_col(
            fill = "grey40",
            width = 0.7
        ) +
        ggplot2::theme_classic(base_size = 11) +
        ggplot2::labs(
            title = "A. Absolute codon usage by amino acid",
            y = "Codon frequency",
            x = NULL
        )

    p2 <- ggplot2::ggplot(df_rscu, ggplot2::aes(x = codon, y = rscu)) +
        ggplot2::geom_col(
            fill = "grey30",
            width = 0.8
        ) +
        ggplot2::geom_hline(
            yintercept = 1,
            linetype = "dashed",
            linewidth = 0.3,
            colour = "black"
        ) +
        ggplot2::facet_wrap(
            ~aa,
            scales = "free_x",
            ncol = 5
        ) +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 90, size = 6),
            axis.title = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(face = "bold")
        ) +
        ggplot2::labs(
            title = "B. Relative synonymous codon usage (RSCU)",
            subtitle = "Dashed line indicates equal synonymous usage (RSCU = 1)"
        )

    p3 <- ggplot2::ggplot(df_rscu, ggplot2::aes(x = aa, y = rscu)) +
        ggplot2::geom_boxplot(
            fill = "grey85",
            width = 0.6,
            outlier.shape = NA
        ) +
        ggplot2::geom_jitter(
            width = 0.15,
            size = 1,
            alpha = 0.6
        ) +
        ggplot2::geom_hline(
            yintercept = 1,
            linetype = "dashed",
            linewidth = 0.3
        ) +
        ggplot2::theme_classic(base_size = 11) +
        ggplot2::labs(
            title = "C. Distribution of RSCU across amino acids",
            y = "RSCU",
            x = NULL
        )

    patchwork::wrap_plots(
        p1, p2, p3,
        design = "
    AB
    CB
    "
    )
}
