#' Plot β-diversity densities within vs across cave types
#'
#' @param comm   community matrix/data.frame (samples in rows, taxa in cols)
#' @param cave factor/character vector of cave labels (length = nrow(comm)).
#' @param beta_obj list of distance objects from BAT::beta() (Btotal, Brepl, Brich)
#'                 or betapart::beta.pair() (beta.sor, beta.sim, beta.sne).
#' @param show_medians logical; add median vlines per group.
#' @param adjust numeric; density bandwidth adjust passed to geom_density().
#' @param facet_cols integer; number of columns for facet layout (1 => 1×k, k => k×1).
#' @param include_total logical; include Total β panel.
#' @param include_nestedness logical; include Nestedness panel.
#' @param fill_values named colors for "Within cave" and "Across caves".
#' @param alpha numeric; fill alpha for densities.
#' @param title character; plot title.
#' @param subtitle character; optional subtitle.
#' @return ggplot object (with raw long data attached as attr "data")
beta_density_plot <- function(
    comm, cave, beta_obj = NULL,
    abund = TRUE,                    # used only if beta_obj is NULL and we compute BAT::beta()
    show_medians = TRUE,
    adjust = 1.2,
    facet_cols = 1,
    include_total = TRUE,
    include_nestedness = TRUE,
    fill_values = c("Within habitat" = "brown",
                    "Across habitats" = "orange"),
    alpha = 0.35,
    title = "",
    subtitle = NULL
) {
  # ---- dependencies check (optional, but nicer errors)
  stopifnot(requireNamespace("dplyr", quietly = TRUE))
  stopifnot(requireNamespace("tibble", quietly = TRUE))
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  
  # ---- coerce comm to numeric matrix
  comm <- as.data.frame(comm, check.names = FALSE)
  # drop non-numeric columns if any (common if Sample column is still there)
  is_num <- vapply(comm, is.numeric, logical(1))
  if (!all(is_num)) comm <- comm[, is_num, drop = FALSE]
  comm <- as.matrix(comm)
  
  stopifnot(nrow(comm) > 1, ncol(comm) > 0)
  
  # Ensure sample names
  if (is.null(rownames(comm))) rownames(comm) <- sprintf("S%03d", seq_len(nrow(comm)))
  sample_ids <- rownames(comm)
  
  # Align cave to comm
  if (is.null(names(cave))) {
    if (length(cave) != nrow(comm)) stop("Un-named 'cave' must have length nrow(comm).")
    cave <- factor(cave); names(cave) <- sample_ids
  } else {
    cave <- factor(cave[sample_ids])
  }
  
  melt_dist <- function(d, label) {
    if (inherits(d, "dist")) {
      labs <- attr(d, "Labels"); if (is.null(labs)) labs <- sample_ids
      m <- as.matrix(d); rownames(m) <- colnames(m) <- labs
    } else {
      m <- as.matrix(d)
      if (is.null(rownames(m)) || is.null(colnames(m))) {
        rownames(m) <- colnames(m) <- sample_ids
      }
    }
    if (!all(sample_ids %in% rownames(m)) || !all(sample_ids %in% colnames(m))) {
      stop("Distance component '", label, "' does not contain all samples from 'comm'.")
    }
    m <- m[sample_ids, sample_ids, drop = FALSE]
    idx <- which(upper.tri(m), arr.ind = TRUE)
    tibble::tibble(
      sample1 = rownames(m)[idx[,1]],
      sample2 = colnames(m)[idx[,2]],
      value   = as.numeric(m[idx]),
      metric  = label
    )
  }
  
  pieces_all <- list()
  
  # ---- If beta_obj is NULL: compute BAT::beta
  if (is.null(beta_obj)) {
    if (!requireNamespace("BAT", quietly = TRUE)) {
      stop("beta_obj is NULL, but package 'BAT' is not installed/available.")
    }
    beta_obj <- BAT::beta(comm, abund = abund)
  }
  
  # ---- If beta_obj is a dist (e.g. Bray): plot just total
  if (inherits(beta_obj, "dist")) {
    if (!include_total) stop("beta_obj is a 'dist' but include_total = FALSE -> nothing to plot.")
    pieces_all[["Total β"]] <- beta_obj
  } else {
    # ---- Otherwise: expect BAT::beta() or betapart::beta.pair()
    have <- names(beta_obj)
    
    if (any(c("Btotal","Brepl","Brich") %in% have)) {
      if (include_total && "Btotal" %in% have) pieces_all[["Total β"]] <- beta_obj$Btotal
      if ("Brepl" %in% have && "Btotal" %in% have) pieces_all[["Turnover"]] <- beta_obj$Brepl / beta_obj$Btotal
      if (include_nestedness && "Brich" %in% have && "Btotal" %in% have) pieces_all[["Nestedness"]] <- beta_obj$Brich / beta_obj$Btotal
      
    } else if (any(c("beta.sor","beta.sim","beta.sne") %in% have)) {
      if (include_total && "beta.sor" %in% have) pieces_all[["Total β (Sørensen)"]] <- beta_obj$beta.sor
      if ("beta.sim" %in% have && "beta.sor" %in% have) pieces_all[["Turnover (Simpson)"]] <- beta_obj$beta.sim / beta_obj$beta.sor
      if (include_nestedness && "beta.sne" %in% have && "beta.sor" %in% have) pieces_all[["Nestedness"]] <- beta_obj$beta.sne / beta_obj$beta.sor
      
    } else {
      stop("Unrecognized 'beta_obj'. Pass a 'dist', BAT::beta() output, or betapart::beta.pair() output.")
    }
  }
  
  if (length(pieces_all) == 0) stop("No metrics to plot. Check include_* flags and beta_obj contents.")
  
  beta_long <- dplyr::bind_rows(lapply(names(pieces_all), function(lbl) {
    melt_dist(pieces_all[[lbl]], lbl)
  })) |>
    dplyr::mutate(
      h1 = cave[sample1],
      h2 = cave[sample2],
      pair_type = dplyr::if_else(h1 == h2, "Within habitat", "Across habitats"),
      metric = factor(metric, levels = c(
        grep("^Total", names(pieces_all), value = TRUE),
        grep("Turnover", names(pieces_all), value = TRUE),
        grep("Nestedness", names(pieces_all), value = TRUE)
      ))
    )
  
  meds <- beta_long |>
    dplyr::group_by(metric, pair_type) |>
    dplyr::summarise(median = stats::median(value, na.rm = TRUE), .groups = "drop")
  
  p <- ggplot2::ggplot(beta_long, ggplot2::aes(x = value, fill = pair_type)) +
    ggplot2::geom_density(alpha = alpha, adjust = adjust, na.rm = TRUE, color = NA) +
    { if (show_medians)
      ggplot2::geom_vline(
        data = meds,
        ggplot2::aes(xintercept = median, linetype = pair_type),
        linewidth = 0.4, show.legend = TRUE
      )
    } +
    ggplot2::facet_wrap(~ metric, ncol = facet_cols, scales = "free_y") +
    ggplot2::scale_fill_manual(values = fill_values, drop = FALSE) +
    ggplot2::scale_linetype_manual(values = c("Within habitat" = "dashed",
                                              "Across habitats" = "solid")) +
    ggplot2::labs(
      x = "Pairwise β-diversity", y = "Density",
      fill = "", linetype = "",
      title = title, subtitle = subtitle
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      legend.position = "top",
      legend.text = element_text(size = 14),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = element_text(size = 12)
    )
}
