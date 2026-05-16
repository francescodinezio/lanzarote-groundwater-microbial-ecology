plot_nmds <- function(comm_matrix, stations, palette, plot_title,
                      distance = "jaccard",
                      point_size = 1.2,
                      label_size = 2.3,
                      label_vjust = -0.5) {
  
  # Reorder rows to match station metadata
  comm_matrix <- comm_matrix[match(rownames(stations), rownames(comm_matrix)), , drop = FALSE]
  stopifnot(all(rownames(comm_matrix) == rownames(stations)))
  
  # Colors by Type
  stations$Type <- factor(stations$Type, levels = names(palette))
  
  # Run NMDS
  nmds <- vegan::metaMDS(comm_matrix, distance = distance, trymax = 100, autotransform = FALSE)
  scores <- as.data.frame(vegan::scores(nmds, display = "sites"))
  scores$Sample <- rownames(scores)
  scores$Type <- stations$Type
  
  stress_txt <- paste0("Stress = ", sprintf("%.3f", nmds$stress))
  
  # Plot
  ggplot2::ggplot(scores, ggplot2::aes(x = NMDS1, y = NMDS2, color = Type)) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::geom_text(ggplot2::aes(label = Sample), vjust = label_vjust, size = label_size) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::labs(title = plot_title) +
    ggplot2::annotate("text", x = Inf, y = Inf, label = stress_txt,
                      hjust = 1.2, vjust = 3.5, size = 5.5) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14)
    )
}