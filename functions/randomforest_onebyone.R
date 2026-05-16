# deps used below
library(randomForest)
library(dplyr)
library(rlang)

# Optional for AUC (set to FALSE if you don't want pROC)
USE_PRC_AUC <- TRUE
if (USE_PRC_AUC) {
  if (!requireNamespace("pROC", quietly = TRUE)) USE_PRC_AUC <- FALSE
}

# Run RF for each habitat vs the rest, for each taxonomic rank
run_rf_one_vs_rest_by_rank <- function(comm.tax.list, stations,
                                       target_col = "Type.1",
                                       ranks = c("Order","Family","Genus","Species"),
                                       ntree = 5000,
                                       balance = TRUE,       # class balancing
                                       plot = TRUE,          # varImp plots
                                       seed = 123) {

  stopifnot(target_col %in% names(stations))
  set.seed(seed)

  habitats_all <- levels(factor(stations[[target_col]]))

  models      <- list()   # models[[rank]][[habitat]]
  importance  <- list()   # importance[[rank]][[habitat]] (data.frame)
  metrics     <- list()   # metrics[[rank]] -> data.frame with habitat, OOB, AUC (if available)

  for (rank in ranks) {
    rf.data <- comm.tax.list[[rank]]
    if (is.null(rf.data) || ncol(rf.data) == 0) {
      message("Skipping rank ", rank, ": empty matrix.")
      next
    }

    # Align rows to stations (safety)
    rf.data <- rf.data[match(rownames(stations), rownames(rf.data)), , drop = FALSE]
    stopifnot(all(rownames(rf.data) == rownames(stations)))

    # Clean column names
    names(rf.data) <- make.names(names(rf.data))

    rank_models     <- list()
    rank_importance <- list()
    rank_metrics    <- list()

    for (hab in habitats_all) {
      # Binary response: focal vs other
      y <- factor(ifelse(stations[[target_col]] == hab, hab, "other"), levels = c("other", hab))
      dat <- data.frame(y = y, rf.data, check.names = FALSE)

      # Skip if only one class present (edge case)
      if (length(unique(y)) < 2) {
        message("Skipping ", rank, " – habitat ", hab, ": only one class present.")
        next
      }

      # Balanced sampling if requested
      sampsize <- NULL
      strata   <- NULL
      if (balance) {
        n_other <- sum(y == "other")
        n_hab   <- sum(y == hab)
        n_min   <- min(n_other, n_hab)
        sampsize <- c(other = n_min, setNames(n_min, hab))
        strata   <- y
      }

      # Fit model
      rf_fit <- randomForest(
        y ~ .,
        data       = dat,
        ntree      = ntree,
        importance = TRUE,
        strata     = strata,
        sampsize   = sampsize
      )

      # Variable importance (sorted)
      imp <- as.data.frame(importance(rf_fit))
      imp$taxon <- rownames(imp)
      imp <- imp[order(imp$MeanDecreaseAccuracy, decreasing = TRUE), ]
      rownames(imp) <- NULL

      # OOB error
      oob <- rf_fit$err.rate[ntree, "OOB"]

      # AUC (optional; uses OOB-style preds via inbag votes)
      auc <- NA_real_
      if (USE_PRC_AUC && ncol(rf_fit$votes) == 2) {
        # probability of being the focal habitat
        p_hat <- rf_fit$votes[, colnames(rf_fit$votes) == hab]
        try({
          auc <- pROC::roc(response = y, predictor = p_hat, quiet = TRUE)$auc[[1]]
        }, silent = TRUE)
      }

      # Store
      rank_models[[hab]]     <- rf_fit
      rank_importance[[hab]] <- imp
      rank_metrics[[length(rank_metrics) + 1]] <-
        data.frame(rank = rank, habitat = hab, n_focal = sum(y == hab),
                   n_other = sum(y == "other"), OOB = oob, AUC = auc, row.names = NULL)

      # Optional plot
      if (plot) {
        main_txt <- paste0("Top 20 taxa – ", rank, " | ", hab, " vs rest")
        try(varImpPlot(rf_fit, n.var = 20, type = 1, main = main_txt), silent = TRUE)
      }
    }

    models[[rank]]     <- rank_models
    importance[[rank]] <- rank_importance
    metrics[[rank]]    <- do.call(rbind, rank_metrics)
  }

  list(models = models, importance = importance, metrics = metrics)
}



bind_metrics <- function(rf_obj, group) {
  # If rf_obj or metrics missing, return empty tibble
  if (is.null(rf_obj) || is.null(rf_obj$metrics)) return(tibble())
  
  mets <- rf_obj$metrics
  
  # Keep only non-null, non-empty elements
  keep_idx <- !vapply(mets, is.null, logical(1)) &
    vapply(mets, function(x) !is.null(x) && NROW(x) > 0, logical(1))
  mets <- mets[keep_idx]
  
  if (length(mets) == 0) return(tibble())
  
  # Some implementations store one data.frame per rank,
  # others store a list of data.frames per rank. Normalize both cases.
  norm <- imap(mets, function(df_or_list, rank_name) {
    if (is.data.frame(df_or_list)) {
      out <- df_or_list
    } else if (is.list(df_or_list)) {
      out <- bind_rows(df_or_list)  # flatten list of dfs for this rank
    } else {
      return(NULL)
    }
    # Ensure required columns exist
    req <- c("rank","habitat","n_focal","n_other","OOB","AUC")
    for (nm in req) if (!nm %in% names(out) && nm == "rank") out$rank <- rank_name
    if (!"rank" %in% names(out)) out$rank <- rank_name
    out
  })
  
  norm <- compact(norm)                # drop NULLs
  if (length(norm) == 0) return(tibble())
  df <- bind_rows(norm)
  df$group <- group
  df
}
