dataframe2fas <- function (x, file = NULL) {
    
    # Check if file name is provided
    if (is.null(file)) {
        "File name to save the results have to be given."
    }

    # Ensure input is a dataframe
    if (!is.data.frame(x)) {
        x <- as.data.frame(x)
    }

    # Handle single-column dataframes
    if (ncol(x) == 1) {
        x <- cbind(rownames(x), x)
        colnames(x) <- c("Seqnames", "squence")
    }

    # Validate correct number of columns
    if (!ncol(x) == 2) {
        stop("Wrong dimention: the input dataframe must be \n either in 1 or 2 dimentions.")
    }

    # Convert columns to character
    dnanames <- as.character(x[, 1])
    dnas <- as.character(x[, 2])

    # Construct FASTA format
    fasta_content <- c()
    for (i in seq_along(dnanames)) {
        fasta_content <- c(fasta_content, paste0(">", dnanames[i]), dnas[i])
    }
    
    # Write to file
    writeLines(fasta_content, file)
    
    # Convert to FASTA format using as.fasta (if applicable)
    if (exists("as.fasta")) {
        return(as.fasta(fasta_content))
    } else {
        return(fasta_content)
    }
}