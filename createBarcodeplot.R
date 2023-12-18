createBarcodeRankPlot <- function(bcrank, all_cells, trim_cells, save_path) {
  # Filter cells
  cells <- all_cells %in% trim_cells
  to_keep <- !duplicated(bcrank$total)
  
  # Create data for graph
  graph_data <- data.frame(Rx = bcrank$rank, Tx = bcrank$total, cells = cells)
  graph_data <- graph_data[to_keep, ]
  
  # Subset the data
  subset_df <- subset(graph_data, graph_data$Tx > 0)

  # Save the plot
  png(save_path, width = 700, height = 500)
  
  # Create a scatter plot
  plot(subset_df$Rx, subset_df$Tx,
       col = ifelse(subset_df$cell, "#0000FF", "#ADD8E6"),
       cex = 1, pch = 21,
       log = "xy", xlab = "Barcodes", ylab = "UMI counts",
       xaxt = "n", yaxt = "n")
  
  # Add horizontal lines and annotations
  abline(h = bcrank@metadata$knee, lty = 2, col = '#FF00FF', lwd = 2)
  text(max(subset_df$Rx) * 0.6, bcrank@metadata$knee + 10000, "Knee", col = "#FF0090", cex = 1)
  
  abline(h = bcrank@metadata$inflection, lty = 2, col = "#FF9966", lwd = 2)
  text(max(subset_df$Rx) * 0.6 , stats@metadata$inflection + 500, "Inflection", col = "#FF8040", cex = 1)
  
  # Set log-scale axis labels
  axis(1, at = 10^seq(floor(log10(min(subset_df$Rx))),
                      ceiling(log10(max(subset_df$Rx)))), 
       labels = prettyNum(10^seq(floor(log10(min(subset_df$Rx))), 
                                 ceiling(log10(max(subset_df$Rx)))), big.mark = ",", scientific = FALSE))
  
  axis(2, at = 10^seq(floor(log10(min(subset_df$Tx))),
                      ceiling(log10(max(subset_df$Tx)))), 
       labels = prettyNum(10^seq(floor(log10(min(subset_df$Tx))), 
                                 ceiling(log10(max(subset_df$Tx)))), big.mark = ",", scientific = FALSE))
  
  # Add legend and title
  legend("bottomleft", c("Cells", "Background"),
         col = c("blue", "lightblue"), pch = c(16,16), cex = 1.2)
  title("Barcode Rank Plot")

  dev.off()
}

# Example usage
#createBarcodeRankPlot(bcrank, all_cells, trim_cells, "path/to/save/plot.png")
