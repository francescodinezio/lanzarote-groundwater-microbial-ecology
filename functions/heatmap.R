mk_heat <- function(mod, title = "",
                    ord = c("cave","pool","salt","sea","well","pond"),
                    disp_lab = c(cave="cave", pool="anc.p", salt="sltw", sea="sea", well="well", pond="spr"),
                    type_var = "Type",
                    digits = 2,
                    alpha_sig = 1,
                    text_white_threshold = 0.06,
                    label_size = 6,
                    axis.text.size = 20,
                    legend.text.size = 20,
                    legend.title.size = 24,
                    limit = 0.3,
                    tile_ratio = NULL,
                    show_x = TRUE,
                    show_y = TRUE) {
  
  # ---- packages ------------------------------------------------------------
  # richiede: emmeans, dplyr, tidyr, ggplot2, scales, grid
  
  # ---- checks --------------------------------------------------------------
  if (!all(ord %in% names(disp_lab))) {
    stop("disp_lab must have names for all ord levels. Missing: ",
         paste(setdiff(ord, names(disp_lab)), collapse = ", "))
  }
  
  term_formula <- stats::as.formula(paste0("~ ", type_var))
  
  # ---- emmeans -------------------------------------------------------------
  emm <- try(emmeans::emmeans(mod, term_formula), silent = TRUE)
  
  if (inherits(emm, "try-error")) {
    dn <- try(as.character(stats::getCall(mod)$data), silent = TRUE)
    dat <- if (!inherits(dn, "try-error") && exists(dn, envir = parent.frame())) {
      get(dn, envir = parent.frame())
    } else {
      NULL
    }
    
    if (!is.null(dat)) {
      rg <- try(emmeans::ref_grid(mod, data = dat), silent = TRUE)
      if (!inherits(rg, "try-error")) {
        emm <- try(emmeans::emmeans(rg, term_formula), silent = TRUE)
      }
    }
    
    if (inherits(emm, "try-error")) {
      stop(
        "emmeans failed for '", type_var, "'. Error: ",
        tryCatch(attr(emm, "condition")$message, error = function(e) "unknown")
      )
    }
  }
  
  # ---- response-scale means ------------------------------------------------
  pr <- summary(emm, type = "response") |> as.data.frame()
  
  prob_col <- intersect(c("prob", "p", "response", "rate", "emmean"), names(pr))
  if (length(prob_col) == 0) {
    stop("No response-scale probability column found in emmeans summary.")
  }
  
  probs <- pr |>
    dplyr::rename(Type = !!type_var) |>
    dplyr::mutate(Type = trimws(as.character(Type))) |>
    dplyr::filter(!is.na(Type)) |>
    dplyr::distinct(Type, .keep_all = TRUE) |>
    dplyr::select(Type, prob = dplyr::all_of(prob_col[1]))
  
  levs <- unique(probs$Type)
  if (length(levs) < 2) {
    stop("Not enough levels to plot. Levels found: ", paste(levs, collapse = ", "))
  }
  
  # ---- pairwise contrasts (make symmetric) ---------------------------------
  pair_df <- as.data.frame(
    summary(emmeans::contrast(emm, method = "pairwise"), adjust = "tukey")
  )
  
  if (!"contrast" %in% names(pair_df)) {
    stop("Internal: expected 'contrast' column not found.")
  }
  
  ct_raw <- pair_df |>
    tidyr::separate(contrast, into = c("i", "j"), sep = " - ", remove = TRUE) |>
    dplyr::mutate(
      i = trimws(i),
      j = trimws(j)
    )
  
  # original orientation
  ct1 <- ct_raw |>
    dplyr::left_join(probs, by = c("i" = "Type")) |>
    dplyr::rename(p_i = prob) |>
    dplyr::left_join(probs, by = c("j" = "Type")) |>
    dplyr::rename(p_j = prob) |>
    dplyr::mutate(
      diff_prob = p_i - p_j,
      sig = !is.na(p.value) & p.value < 0.05
    )
  
  # mirrored orientation
  ct2 <- ct_raw |>
    dplyr::transmute(i_sw = j, j_sw = i, estimate, SE, df, z.ratio, p.value) |>
    dplyr::rename(i = i_sw, j = j_sw) |>
    dplyr::left_join(probs, by = c("i" = "Type")) |>
    dplyr::rename(p_i = prob) |>
    dplyr::left_join(probs, by = c("j" = "Type")) |>
    dplyr::rename(p_j = prob) |>
    dplyr::mutate(
      diff_prob = p_i - p_j,
      sig = !is.na(p.value) & p.value < 0.05
    )
  
  ct_sym <- dplyr::bind_rows(ct1, ct2) |>
    dplyr::filter(i != j)
  
  # ---- complete off-diagonal grid ------------------------------------------
  ord_present <- ord[ord %in% levs]
  
  grid <- tidyr::expand_grid(i = ord_present, j = ord_present) |>
    dplyr::filter(i != j)
  
  ct <- grid |>
    dplyr::left_join(ct_sym, by = c("i", "j"))
  
  # ---- labels & aesthetics -------------------------------------------------
  fmt <- paste0("%.", digits, "f")
  
  ct <- ct |>
    dplyr::mutate(
      label = dplyr::if_else(!is.na(sig) & sig, sprintf(fmt, diff_prob), "-"),
      text_col_key = dplyr::case_when(
        is.na(sig) ~ "nsg",
        !sig ~ "nsg",
        abs(diff_prob) > text_white_threshold ~ "white",
        TRUE ~ "black"
      )
    )
  
  ct$i <- factor(ct$i, levels = ord_present)
  ct$j <- factor(ct$j, levels = ord_present)
  
  # ---- dynamic axis theme --------------------------------------------------
  theme_axes <- ggplot2::theme(
    axis.text.x = if (show_x) {
      ggplot2::element_text(size = axis.text.size, angle = 45, hjust = 1, vjust = 1)
    } else {
      ggplot2::element_blank()
    },
    axis.text.y = if (show_y) {
      ggplot2::element_text(size = axis.text.size)
    } else {
      ggplot2::element_blank()
    },
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.line.x  = ggplot2::element_blank(),
    axis.line.y  = ggplot2::element_blank()
  )
  
  # ---- plot ----------------------------------------------------------------
  p <- ggplot2::ggplot(ct, ggplot2::aes(i, j, fill = diff_prob)) +
    ggplot2::geom_tile(color = "grey90", alpha = alpha_sig) +
    ggplot2::geom_text(
      ggplot2::aes(label = label, color = text_col_key),
      fontface = "bold",
      size = label_size,
      na.rm = TRUE
    ) +
    ggplot2::scale_fill_gradient2(
      name = "\u0394 prob.",
      limits = c(-limit, limit),
      midpoint = 0,
      breaks = c(-limit, 0, limit),
      labels = scales::number_format(accuracy = 0.01),
      oob = scales::squish
    ) +
    ggplot2::scale_color_manual(
      values = c(white = "white", black = "black", nsg = "grey30"),
      guide = "none"
    ) +
    ggplot2::scale_x_discrete(
      limits = ord_present,
      labels = disp_lab[ord_present],
      drop = FALSE
    ) +
    ggplot2::scale_y_discrete(
      limits = ord_present,
      labels = disp_lab[ord_present],
      drop = FALSE
    ) +
    ggplot2::labs(x = "", y = "", title = title) +
    ggplot2::theme_classic(base_size = 10) +
    theme_axes +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(size = legend.title.size, face = "bold"),
      legend.text  = ggplot2::element_text(size = legend.text.size),
      legend.key.height = grid::unit(0.7, "cm"),
      legend.key.width  = grid::unit(0.7, "cm"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    )
  
  if (!is.null(tile_ratio)) {
    p <- p + ggplot2::coord_fixed(ratio = tile_ratio)
  }
  
  return(p)
}