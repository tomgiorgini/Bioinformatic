getEnrichmentPlot <- function(annotation, type, top_term, thr_pval, dirEnrich){
  
  # Filtra per p-value
  annotation <- annotation[annotation$Adjusted.P.value < thr_pval, ]
  annotation <- annotation[order(annotation$Adjusted.P.value), ]
  
  # Se non ci sono termini significativi, esci silenziosamente
  if(nrow(annotation) == 0){
    message(sprintf("getEnrichmentPlot: nessun termine con Adjusted.P.value < %s, salto plotting.", thr_pval))
    return(invisible(NULL))
  }
  
  # Calcola Gene_count e Gene_ratio
  annotation$Gene_count <- sapply(annotation$Genes, function(x){
    length(strsplit(x, ";")[[1]])
  })
  
  annotation$Gene_ratio <- sapply(annotation$Overlap, function(x){
    parts <- as.numeric(strsplit(x, "/")[[1]])
    parts[1] / parts[2]
  })
  
  # Seleziona i top_term se richiesto
  if(length(top_term) != 0 && top_term <= nrow(annotation)){
    annotation_top <- annotation[1:top_term, ]
  } else {
    annotation_top <- annotation
  }
  
  # Se annotation_top è vuoto, skip plotting
  if(nrow(annotation_top) == 0){
    message("getEnrichmentPlot: dopo selezione top_term, nulla da plottare.")
  } else {
    # Bar plot
    g1 <- ggplot(annotation_top, 
                 aes_string(x = "Gene_count",
                            y = "fct_reorder(Term, Gene_count)",
                            fill = "Adjusted.P.value")) +
      geom_bar(stat="identity") +
      scale_fill_continuous(low="red", high="blue", name="Adjusted.P.value",
                            guide=guide_colorbar(reverse=TRUE)) +
      scale_y_discrete(labels = function(x) str_wrap(x, width=40)) +
      theme_bw(base_size=10) +
      ylab(NULL)
    
    pdf(file.path(dirEnrich, paste0(type, "_barplot.pdf")))
    print(g1)
    dev.off()
    
    # Dot plot
    g2 <- ggplot(annotation_top, 
                 aes_string(x = "Gene_count",
                            y = "fct_reorder(Term, Gene_count)")) +
      geom_point(aes(size = Gene_ratio, color = Adjusted.P.value)) +
      scale_colour_gradient(limits = c(min(annotation_top$Adjusted.P.value),
                                       max(annotation_top$Adjusted.P.value)),
                            low="red", high="blue") +
      theme_bw(base_size=10) +
      scale_y_discrete(labels = function(x) str_wrap(x, width=40)) +
      ylab(NULL)
    
    pdf(file.path(dirEnrich, paste0(type, "_dotplot.pdf")))
    print(g2)
    dev.off()
  }
  
  # Esporta tabella dei risultati (anche se vuota, non dà errore)
  out_file <- file.path(dirEnrich,
                        paste0(type, "_adj_pval_", thr_pval, ".txt"))
  write.table(annotation[, c("Term","Overlap","P.value","Adjusted.P.value",
                             "Gene_count","Gene_ratio","Genes")],
              out_file, sep="\t", quote=FALSE, row.names=FALSE)
}