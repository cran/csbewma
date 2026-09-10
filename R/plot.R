utils::globalVariables(c("r", "UCL", "LCL", "stream", "p_value", "color_status"))
#' Plot CSB-EWMA Control Chart
#'
#' Creates a professional control chart showing EWMA statistic and limits.
#'
#' @param x csb_ewma object from csb_ewma() function
#' @param title Plot title (default = "CSB-EWMA Control Chart")
#' @param show_signal Whether to highlight signal point (default = TRUE)
#' @param ... Additional arguments passed to ggplot
#' @return A ggplot object (invisibly) and displays the plot
#' @export
#' @examples
#' # See csb_ewma() for examples
plot.csb_ewma <- function(x, title = "CSB-EWMA Control Chart",
                          show_signal = TRUE, ...) {
  
  # Check that x is a valid csb_ewma object
  if (!inherits(x, "csb_ewma")) {
    stop("Object must be of class 'csb_ewma'")
  }
  
  # Extract data from the csb_ewma object
  t_values <- 1:x$T_sig
  
  # Create data frame for ggplot
  plot_data <- data.frame(
    t = t_values,
    r = x$r_history,
    UCL = x$UCL_history,
    LCL = x$LCL_history
  )
  
  # Create the base ggplot object with explicit ggplot2:: prefix
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = t)) +
    
    ggplot2::geom_line(ggplot2::aes(y = r, color = "EWMA Statistic"), linewidth = 1) +
    ggplot2::geom_line(ggplot2::aes(y = UCL, color = "UCL"), linetype = "dashed", linewidth = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = LCL, color = "LCL"), linetype = "dashed", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "gray50", linewidth = 0.5) +
    
    ggplot2::labs(
      title = title,
      x = "Time Point (t)",
      y = expression(r[t] ~ "(EWMA Statistic)"),
      color = "Series"
    ) +
    
    ggplot2::scale_color_manual(
      values = c("EWMA Statistic" = "#1F7AB5",
                 "UCL" = "#D62828",
                 "LCL" = "#D62828")
    ) +
    
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold", size = 10),
      legend.text = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = ggplot2::element_text(face = "bold", size = 11),
      axis.text = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "gray90", linewidth = 0.3),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
  
  # Add signal point annotation if a signal was detected
  if (show_signal && x$signal_detected) {

    signal_times_to_plot <- if (!is.null(x$signal_times)) x$signal_times else x$signal_time

    signal_data <- data.frame(
      t = signal_times_to_plot,
      r = x$r_history[signal_times_to_plot]
    )

    p <- p + ggplot2::geom_point(data = signal_data, ggplot2::aes(x = t, y = r),
                                  color = "#D62828", size = 4.5, shape = 18) +

      ggplot2::annotate("text",
                         x = signal_times_to_plot,
                         y = x$r_history[signal_times_to_plot] + 0.15,
                         label = paste("Signal at t =", signal_times_to_plot),
                         color = "#D62828",
                         size = 3.5,
                         fontface = "bold",
                         hjust = 1)
  }
  
  # Display the plot
  print(p)
  
  # Return the plot object invisibly
  invisible(p)
}


#' Direct Plot for CSB-EWMA Results
#'
#' Creates a professional control chart showing EWMA statistic and limits.
#' This function can be called directly without S3 dispatch.
#'
#' @param result csb_ewma object from csb_ewma() or run_csb_ewma()
#' @param title Plot title (default = "CSB-EWMA Control Chart")
#' @param show_signal Whether to highlight signal point (default = TRUE)
#' @return A ggplot object (invisibly) and displays the plot
#' @export
plot_csb_ewma_direct <- function(result, title = "CSB-EWMA Control Chart",
                                 show_signal = TRUE) {
  
  # Extract data from the result object
  t_values <- 1:result$T_sig
  
  # Create data frame for plotting
  plot_data <- data.frame(
    t = t_values,
    r = result$r_history,
    UCL = result$UCL_history,
    LCL = result$LCL_history
  )
  
  # Create the base ggplot object with explicit ggplot2:: prefix
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = t)) +
    
    ggplot2::geom_line(ggplot2::aes(y = r, color = "EWMA Statistic"), linewidth = 1) +
    ggplot2::geom_line(ggplot2::aes(y = UCL, color = "UCL"), linetype = "dashed", linewidth = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = LCL, color = "LCL"), linetype = "dashed", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "gray50", linewidth = 0.5) +
    
    ggplot2::labs(
      title = title,
      x = "Time Point (t)",
      y = expression(r[t] ~ "(EWMA Statistic)"),
      color = "Series"
    ) +
    
    ggplot2::scale_color_manual(
      values = c("EWMA Statistic" = "#1F7AB5",
                 "UCL" = "#D62828",
                 "LCL" = "#D62828")
    ) +
    
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold", size = 10),
      legend.text = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = ggplot2::element_text(face = "bold", size = 11),
      axis.text = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "gray90", linewidth = 0.3),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
  
  # Add signal point if detected
  if (show_signal && result$signal_detected) {

    signal_times_to_plot <- if (!is.null(result$signal_times)) result$signal_times else result$signal_time

    p <- p + ggplot2::geom_point(data = data.frame(t = signal_times_to_plot,
                                                     r = result$r_history[signal_times_to_plot]),
                                  ggplot2::aes(x = t, y = r),
                                  color = "#D62828", size = 4.5, shape = 18) +

      ggplot2::annotate("text",
                         x = signal_times_to_plot,
                         y = result$r_history[signal_times_to_plot] + 0.15,
                         label = paste("Signal at t =", signal_times_to_plot),
                         color = "#D62828",
                         size = 3.5,
                         fontface = "bold",
                         hjust = 1)
  }
  
  # Display the plot
  print(p)
  
  # Return the plot object invisibly
  invisible(p)
}


#' Plot Flagged Streams Bar Chart
#'
#' Creates a bar plot showing -log10(p-values) for each stream.
#' Flagged streams appear in red, others in gray.
#'
#' @param flagged_results Output data frame from identify_ooc() function
#' @param alpha Significance level for reference line (default = 0.05)
#' @param title Plot title (default = "Stream P-values")
#' @return A ggplot object (invisibly) and displays the plot
#' @export
plot_flagged_streams <- function(flagged_results, alpha = 0.05,
                                 title = "Stream P-values") {
  
  # Create a new column for bar colors based on flagged status
  flagged_results$color_status <- ifelse(flagged_results$flagged,
                                         "Flagged",
                                         "Not Flagged")
  
  # Create the bar plot with explicit ggplot2:: prefix
  p <- ggplot2::ggplot(flagged_results, ggplot2::aes(x = factor(stream),
                                                     y = -log10(p_value),
                                                     fill = color_status)) +
    
    ggplot2::geom_bar(stat = "identity", width = 0.7) +
    
    ggplot2::geom_hline(yintercept = -log10(alpha),
                        linetype = "dashed",
                        color = "#D62828",
                        linewidth = 0.8) +
    
    ggplot2::annotate("text",
                      x = nrow(flagged_results),
                      y = -log10(alpha) + 0.1,
                      label = paste("alpha =", alpha),
                      color = "#D62828",
                      size = 3,
                      hjust = 1) +
    
    ggplot2::labs(
      title = title,
      x = "Stream Number",
      y = expression(-log[10](p-value)),
      fill = "Status"
    ) +
    
    ggplot2::scale_fill_manual(values = c("Flagged" = "#D62828",
                                          "Not Flagged" = "#A9A9A9")) +
    
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold", size = 10),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = ggplot2::element_text(face = "bold", size = 11),
      axis.text = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
  
  # Display the plot
  print(p)
  
  # Return the plot object invisibly
  invisible(p)
}


#' Combined Diagnostic Dashboard
#'
#' Creates a combined plot showing both the CSB-EWMA control chart
#' and the flagged streams bar plot side by side or stacked.
#'
#' @param chart_result csb_ewma object from csb_ewma() function
#' @param flagged_results Output from identify_ooc() function
#' @param layout Either "side" for side-by-side or "stacked" for vertical
#' @return A combined ggplot object (invisibly) and displays the plot
#' @export
plot_chart_with_flagged <- function(chart_result, flagged_results,
                                    layout = "side") {
  
  # Create the control chart plot
  p1 <- plot_csb_ewma_direct(chart_result, show_signal = TRUE)
  
  # Create the flagged streams bar plot
  p2 <- plot_flagged_streams(flagged_results, alpha = 0.05)
  
  # Combine plots based on layout choice using patchwork
  if (layout == "side") {
    combined_plot <- p1 + p2 + patchwork::plot_layout(ncol = 2, widths = c(1, 0.8))
  } else {
    combined_plot <- p1 / p2 + patchwork::plot_layout(nrow = 2, heights = c(1, 0.8))
  }
  
  # Add a common title
  combined_plot <- combined_plot +
    patchwork::plot_annotation(
      title = "CSB-EWMA Diagnostic Dashboard",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16))
    )
  
  # Display the combined plot
  print(combined_plot)
  
  # Return invisibly
  invisible(combined_plot)
}
