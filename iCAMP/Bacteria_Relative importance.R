library(ggplot2)
library(cowplot)

# 
df <- read.delim("Bacteria_Relative importance.txt",
                 header = TRUE,
                 sep = "\t",
                 check.names = FALSE)

# 
processes <- c("Deterministic", "Stochastic", "HeS", "HoS", "DL", "DR", "HD")

dir.create("Scatter_process_plots", showWarnings = FALSE)

plot_one_process <- function(proc){
  
  ck_col <- paste0(proc, "_CK")
  tr_col <- paste0(proc, "_TR")
  
  plot_df <- df[, c(ck_col, tr_col)]
  colnames(plot_df) <- c("CK", "TR")
  plot_df <- na.omit(plot_df)
  
  plot_df$Group <- ifelse(plot_df$TR > plot_df$CK,
                          "TR > CK",
                          "CK > TR")
  
  group_counts <- table(plot_df$Group)
  
  n_tr <- ifelse("TR > CK" %in% names(group_counts),
                 as.integer(group_counts["TR > CK"]), 0)
  n_ck <- ifelse("CK > TR" %in% names(group_counts),
                 as.integer(group_counts["CK > TR"]), 0)
  
  annotation_text <- paste0("TR > CK: ", n_tr,
                            "    CK > TR: ", n_ck)
  
  p <- ggplot(plot_df, aes(x = CK, y = TR, color = Group)) +
    geom_point(size = 5, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "grey40") +
    scale_x_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = c(
      "TR > CK" = "#448DCD",
      "CK > TR" = "#F7AF34"
    )) +
    labs(x = paste0(proc, "_CK"),
         y = paste0(proc, "_TR"),
         color = NULL,
         title = proc) +
    annotate("text",
             x = 0.05, y = 0.95,
             label = annotation_text,
             hjust = 0, vjust = 1,
             size = 4.2,
             fontface = "bold") +
    coord_fixed() +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ggsave(filename = paste0("Scatter_process_plots/Bacteria_", proc, "_Scatter.pdf"),
         plot = ggdraw(p),
         width = 6,
         height = 6,
         units = "in")
  
  return(p)
}

plot_list <- lapply(processes, plot_one_process)

