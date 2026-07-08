

################################################################################################################### PCoA


library(ggplot2)

## Read and merge data
metrics_df <- read.delim("guild_abundance.xls.txt",  check.names = FALSE)
metadata <- read.delim("Fungi_metadata_raw.txt", header = TRUE)

# Merge the two files by the first column (usually SampleID)
# (If the first column names differ, specify using by.x="..." and by.y="...")
net_meta <- merge(metadata, metrics_df, by = 1)

## Non-parametric test for two treatments + BH correction
# Automatically extract all network metric names 
# (Assuming the first column in the CSV is SampleID, and the rest are network metrics)
vars <- setdiff(colnames(metrics_df), colnames(metrics_df)[1])

raw_p <- sapply(vars, function(v) {
  # Extract Group and the corresponding metric column
  dat <- net_meta[, c("Group", v)]
  dat <- dat[complete.cases(dat), ]
  # Wilcoxon rank-sum test
  wilcox.test(dat[[v]] ~ dat$Group, exact = FALSE)$p.value
})

# BH multiple testing correction
p_adj <- p.adjust(raw_p, method = "BH")

# Format the statistical results table
stat_res <- data.frame(
  Index = vars,
  P_raw = raw_p,
  P_adj = p_adj,
  stringsAsFactors = FALSE
)
print(stat_res)

## Export statistical results
write.table(stat_res, "Fungi_Subnetwork_Wilcoxon_BH.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## colors
my_fill <- c("CK" = "#EFB44F", "TR" = "#5E97C9")


## Define plotting function
plot_metrics_df <- function(data, yvar, stat_table, ylab = NULL,
                             out_pdf = NULL, out_png = NULL) {
  # Extract adjusted p-value
  p_use <- stat_table$P_adj[stat_table$Index == yvar]
  p_lab <- paste0("Wilcoxon, BH-adjusted p = ", signif(p_use, 3))
  
  # Remove NA values
  dat <- data[, c("Group", yvar)]
  dat <- dat[complete.cases(dat), ]
  colnames(dat)[2] <- "Value"
  
  # Text positioning for the p-value label
  y_max <- max(dat$Value, na.rm = TRUE)
  y_min <- min(dat$Value, na.rm = TRUE)
  y_pos <- y_max + (y_max - y_min) * 0.08
  
  p <- ggplot(dat, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot(width = 0.55, alpha = 0.9, outlier.shape = NA, linewidth = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, fill = "black", color = "black") +
    scale_fill_manual(values = my_fill) +
    labs(x = NULL, y = ifelse(is.null(ylab), yvar, ylab)) +
    annotate("text", x = 1.5, y = y_pos, label = p_lab, size = 5) +
    theme_classic(base_size = 16) +
    theme(legend.position = "none",
          axis.line = element_line(linewidth = 0.5), 
          axis.ticks = element_line(linewidth = 0.5),
          axis.text.x = element_text(angle = 25, hjust = 1), 
          plot.margin = margin(10, 15, 10, 10)) +
    coord_cartesian(ylim = c(y_min, y_pos * 1.02))
  
  print(p)
  if (!is.null(out_pdf)) {
    ggsave(out_pdf, p, width = 6, height = 5)
  }
  return(p)
}


## Batch plotting and exporting
# Use lapply to iterate over all network metrics, automatically plotting and saving
# This handles any number of metrics in the CSV automatically with one run
plot_list <- lapply(vars, function(v) {
  
  # Automatically generate the corresponding PDF filename
  out_name <- paste0("Fungi_Subnetwork_", v, "_boxplot.pdf")
  
  # Call the plotting function
  plot_metrics_df(
    data = net_meta,
    yvar = v,
    stat_table = stat_res,
    ylab = v,           # Directly use the network metric name as the Y-axis label
    out_pdf = out_name
  )
})

###### Plant Pathogen
p_single <- plot_metrics_df(
  data = net_meta,
  yvar = "Plant Pathogen",       
  stat_table = stat_res,
  ylab = "Plant Pathogen"         
)
p_single
p_single_modified <- p_single + coord_cartesian(ylim = c(0, 0.4))
print(p_single_modified)
ggsave("Fungi_Plant Pathogen.pdf", p_single_modified, width = 6, height = 5)








## ALDEx2

library(ALDEx2)

func <- read.delim("guild_abundance_ALDEx2.txt",
                   header = TRUE,
                   sep = "\t",
                   check.names = FALSE)

metadata <- read.delim("Fungi_metadata_raw.txt",
                       header = TRUE,
                       sep = "\t",
                       check.names = FALSE)

rownames(func) <- func$guild_abundance
func <- func[, -1]

func[] <- lapply(func, function(x) as.numeric(as.character(x)))
func[is.na(func)] <- 0


common_samples <- intersect(colnames(func), metadata$SampleID)
func <- func[, common_samples]
metadata <- metadata[match(common_samples, metadata$SampleID), ]


all(colnames(func) == metadata$SampleID)
keep <- metadata$Group %in% c("CK", "TR")
func <- func[, keep]
metadata <- metadata[keep, ]

func[] <- lapply(func, function(x) {
  as.numeric(as.character(x))
})

func <- func[rowSums(func) > 0, ]

## Group
group <- factor(metadata$Group, levels = c("CK", "TR"))
table(group)

## ALDEx2 analysis
group <- as.character(group)
set.seed(1)
aldex_clr <- aldex.clr(func,
                       conds = group,
                       mc.samples = 128,
                       denom = "all",
                       verbose = TRUE)

aldex_res <- aldex.ttest(aldex_clr,
                         paired.test = FALSE)
aldex_eff <- aldex.effect(aldex_clr)
res <- cbind(aldex_res, aldex_eff)

## Result
res$Function <- rownames(res)

res <- res[, c("Function",
               "we.ep", "we.eBH",
               "wi.ep", "wi.eBH",
               "effect",
               "diff.btw",
               "diff.win")]

res <- res[order(res$we.eBH), ]

res$Significance <- ifelse(res$wi.eBH < 0.001, "***",
                           ifelse(res$wi.eBH < 0.01, "**",
                                  ifelse(res$wi.eBH < 0.05, "*", "ns")))
## Output
write.table(res,
            "ALDEx2_CK_vs_TR_results.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

head(res)

############################################################################## Plot
library(ggplot2)

res_plot <- res

res_plot <- res_plot[order(res_plot$effect), ]
res_plot$Function <- factor(res_plot$Function, levels = res_plot$Function)

res_plot$Significance <- ifelse(res_plot$wi.eBH < 0.001, "***",
                                ifelse(res_plot$wi.eBH < 0.01, "**",
                                       ifelse(res_plot$wi.eBH < 0.05, "*", "")))

p <- ggplot(res_plot, aes(x = effect, y = Function)) +
  geom_bar(stat = "identity",
           aes(fill = effect > 0),
           width = 0.7) +
  
  geom_vline(xintercept = 0, color = "black") +
  
  geom_text(aes(label = Significance),
            hjust = ifelse(res_plot$effect > 0, -0.2, 1.2),
            size = 4) +
  
  scale_fill_manual(values = c("TRUE" = "#F7AF34",   # TR ↑
                               "FALSE" = "#448DCD"), # CK ↑
                    guide = "none") +
  
  labs(x = "Effect size (ALDEx2)",
       y = "") +
  
  theme_classic(base_size = 14)
p
ggsave("effect_barplot.pdf",
       plot = p,
       width = 8,
       height = 6,
       device = cairo_pdf)





################################################################################ Mean and SD calculation


library(dplyr)

metadata <- read.delim("Fungi_metadata_raw.txt", stringsAsFactors = FALSE)
KO <- read.delim("guild_abundance_ALDEx2.txt", check.names = FALSE)

rownames(KO) <- KO[,1]
KO <- KO[,-1]

groups <- unique(metadata$Group_Study)

result_list <- list()

for (g in groups) {
  
  samples <- metadata$SampleID[metadata$Group_Study == g]
  
  samples_exist <- intersect(samples, colnames(KO))
 
  sub_data <- KO[, samples_exist, drop = FALSE]
  
  #  mean and sd
  mean_vec <- apply(sub_data, 1, function(x) {
    if (length(x) == 0) return(NA)
    if (length(x) == 1) return(x)
    mean(x, na.rm = TRUE)
  })
  
  sd_vec <- apply(sub_data, 1, function(x) {
    if (length(x) <= 1) return(NA)
    sd(x, na.rm = TRUE)
  })
  
  # Save
  result_list[[g]] <- data.frame(
    Group = g,
    Feature = rownames(KO),
    Mean = mean_vec,
    SD = sd_vec
  )
}

#  Merge
final_result <- do.call(rbind, result_list)

# Output
write.table(final_result, "guild_abundance_mean_sd.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)















































