summarize_by_taxonomic_rank <- function(com.tax, ranks = c("Order", "Family", "Genus", "Species"), min_reads = 20) {
  comm.tax.list <- list()
  
  for (rank in ranks) {
    # Filter out unclassified
    comm.filtered <- com.tax[com.tax[[rank]] != "Unclassified", ]
    
    # Summarize numeric columns by taxon rank
    comm.sum <- comm.filtered %>%
      group_by(.data[[rank]]) %>%
      summarise(across(where(is.numeric), sum), .groups = "drop") %>%
      as.data.frame()
    
    rownames(comm.sum) <- comm.sum[[rank]]
    comm.sum <- comm.sum[, -1]
    comm.sum <- as.data.frame(t(comm.sum))
    
    # Filter taxa with fewer than `min_reads` reads
    taxon_totals <- colSums(comm.sum)
    comm.sum <- comm.sum[, taxon_totals >= min_reads, drop = FALSE]
    
    # Store in list
    comm.tax.list[[rank]] <- comm.sum
  }
  
  return(comm.tax.list)
}
