habitat_heatmap <- function(df,
                            panel_title = "a",
                            n_show = 10,
                            show_unclassified = FALSE,
                            fill_breaks = c(0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.075, 0.1, 0.4),
                            legend_breaks = c(0, 0.1, 0.2, 0.3, 0.4),
                            fill_colours = c("white", "#d9f0d3", "#a6dba0", "#5ab4ac", "#2b8cbe", "#08589e", "darkblue"),
                            text_white_threshold = 0.01) {
  
  df2 <- df %>%
    dplyr::mutate(
      Group = groups_named[as.character(Sample)],
      Genus = as.character(Genus),
      rel_abund = as.numeric(rel_abund)
    ) %>%
    dplyr::filter(!is.na(Group), !is.na(Genus))
  
  if (!show_unclassified) {
    df2 <- df2 %>% dplyr::filter(Genus != "Unclassified")
  }
  
  genus_sample <- df2 %>%
    dplyr::group_by(Sample, Group, Genus) %>%
    dplyr::summarise(rel_abund = sum(rel_abund, na.rm = TRUE), .groups = "drop")
  
  genus_habitat <- genus_sample %>%
    dplyr::group_by(Group, Genus) %>%
    dplyr::summarise(mean_rel_abund = mean(rel_abund, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Group = factor(Group, levels = hab_order)) %>%
    tidyr::complete(Group, Genus, fill = list(mean_rel_abund = 0))
  
  top_genera <- genus_habitat %>%
    dplyr::group_by(Genus) %>%
    dplyr::summarise(
      mean_abund = mean(mean_rel_abund, na.rm = TRUE),
      best_group = as.character(Group[which.max(mean_rel_abund)][1]),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(mean_abund), Genus) %>%
    dplyr::slice_head(n = n_show)
  
  hm <- genus_habitat %>%
    dplyr::filter(Genus %in% top_genera$Genus) %>%
    dplyr::left_join(top_genera[, c("Genus", "best_group", "mean_abund")], by = "Genus")
  
  genus_order <- hm %>%
    dplyr::distinct(Genus, best_group, mean_abund) %>%
    dplyr::mutate(best_group = factor(best_group, levels = hab_order)) %>%
    dplyr::arrange(best_group, dplyr::desc(mean_abund), Genus) %>%
    dplyr::pull(Genus)
  
  hm <- hm %>%
    dplyr::mutate(
      mean_rel_abund = tidyr::replace_na(mean_rel_abund, 0),
      mean_val = mean_rel_abund,
      label = dplyr::if_else(mean_val < 0.001, "<0.001", sprintf("%.3f", mean_val)),
      text_color = dplyr::if_else(mean_val >= text_white_threshold, "white", "black"),
      Genus = factor(Genus, levels = rev(genus_order)),
      Group = factor(Group, levels = hab_order)
    )
  
  fill_max <- max(c(hm$mean_val, fill_breaks), na.rm = TRUE)
  
  ggplot2::ggplot(hm, ggplot2::aes(x = Group, y = Genus, fill = mean_val)) +
    ggplot2::geom_tile(color = "grey90", linewidth = 0.2) +
    ggplot2::geom_text(
      ggplot2::aes(label = label, colour = text_color),
      fontface = "bold",
      size = 2.8
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::coord_fixed() +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_fill_gradientn(
      colours = fill_colours,
      values = scales::rescale(fill_breaks, from = c(0, fill_max)),
      limits = c(0, fill_max),
      breaks = legend_breaks,
      labels = scales::number_format(accuracy = 0.001),
      na.value = "white",
      oob = scales::squish,
      name = "Mean rel.\nabundance"
    ) +
    ggplot2::labs(title = panel_title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, size = 11),
      axis.text.y = ggplot2::element_text(face = "italic", size = 10),
      plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0),
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 10)
    )
}