

################################################################################################################### PCoA
# ---- Packages ----
library(ggplot2)
library(vegan)

# ---- Inputs ----
otu <- read.delim("Bacteria_otutab_rare.txt", row.names = 1, check.names = FALSE)
metadata <- read.delim("Bacteria_metadata_raw.txt", header = TRUE,  stringsAsFactors = FALSE)
metadata <- metadata[metadata$SampleID %in% colnames(otu), ]
metadata <- metadata[match(colnames(otu), metadata$SampleID), ]
stopifnot(all(metadata$SampleID == colnames(otu)))
metadata$Group <- as.factor(metadata$Group)
metadata$StudyID <- as.factor(metadata$StudyID)

# 2. Bray-Curtis
library(parallelDist)
bray_dist <- parDist( t(otu),  method = "bray", threads = 24)

# 3. PCoA
pcoa <- cmdscale(bray_dist, k = 2, eig = TRUE)
pcoa_points <- as.data.frame(pcoa$points)
pcoa_points$SampleID <- rownames(pcoa_points)
colnames(pcoa_points)[1:2] <- c("PCoA1", "PCoA2")

# Merge metadata
df <- merge(pcoa_points, metadata, by = "SampleID")

# 4. Explaination
eig <- pcoa$eig
eig_pos <- eig[eig > 0]
p1 <- eig[1] / sum(eig_pos) * 100
p2 <- eig[2] / sum(eig_pos) * 100

# 5. PERMANOVA
library(parallel)
cl <- makeCluster(24)
adon <- adonis2(bray_dist ~ StudyID * Group,data = metadata, permutations = 999, parallel = cl)
stopCluster(cl)
print(adon)

r2_group <- adon$R2[1]
p_group  <- adon$`Pr(>F)`[1]

# 6. Color
my_colors <- c("CK" = "#F7AF34",
               "TR" = "#448DCD")
# 7. Plot
p <- ggplot(df, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = Group, shape = Group), size = 3, alpha = 0.85) +
  stat_ellipse(aes(group = Group, color = Group),
               level = 0.95, linetype = "dashed", linewidth = 0.7) +
  scale_color_manual(values = my_colors) +
  xlab(sprintf("PCoA1 (%.2f%%)", p1)) +
  ylab(sprintf("PCoA2 (%.2f%%)", p2)) +
  ggtitle(sprintf(
    paste0("PERMANOVA  Group: R2 = %.2f%%, P = %.3f; "),
    r2_group * 100, p_group)) +
  theme_bw() +
  theme(
    axis.text = element_text(colour = "black", size = 9),
    axis.title = element_text(colour = "black", size = 11),
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank()
  )
p

p1 <- p + coord_cartesian(ylim = c(-0.5, 0.35))
print(p1)

# 8. Save
ggsave("Bacteria_PCoA_PERMANOVA.pdf", plot = p1,
       width = 7, height = 6.5, units = "in")

# 9. Output
write.csv(df[, c("SampleID", "PCoA1", "PCoA2", "Group", "StudyID")],
          "Bacteria_PCoA_Coordinates.csv",
          row.names = FALSE)
cat("PCoA coordinates saved to 'Bacteria_PCoA_Coordinates.csv'\n")




#################################################################################################### Community variability
# 1. Library
library(vegan)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)

# Calculate beta dispersion
bd <- betadisper(bray_dist, metadata$Group) # "Group"

# Plot
plot(bd)
boxplot(bd, main = "Beta Dispersion by Group")

anova(bd)
permutest(bd, permutations = 999)

# Distance to centroid
dispersion_result <- data.frame(
  SampleID = names(bd$distances),
  Group = bd$group,
  DistanceToCentroid = bd$distances
)
# Results
head(dispersion_result)
# Output
write.csv(dispersion_result, "Bacteria_beta_dispersion_result.csv", row.names = FALSE)

# data from betadisper 
df_disp <- data.frame(
  SampleID = names(bd$distances),
  Group    = as.character(bd$group),
  Distance = as.numeric(bd$distances)
)

df_disp$Group <- factor(df_disp$Group, levels = c("CK","TR"))
fill_cols <- c("CK"="#F7AF34","TR"="#448DCD")

perm <- permutest(bd, permutations = 999)
p_perm <- tryCatch(
  as.numeric(perm$tab[1, "Pr(>F)"]),
  error = function(e) NA_real_
)
p_label <- ifelse(
  p_perm < 0.001,
  "p < 0.001",
  paste0("p = ", signif(p_perm, 3))
)

custom_ylim <- c(0.6, 0.8)

y_min <- custom_ylim[1]
y_max <- custom_ylim[2]

y_line <- y_max * 0.95
y_text <- y_max * 0.97

# Plot
p <- ggplot(df_disp, aes(x = Group, y = Distance, fill = Group)) +
  geom_boxplot(
    width = 0.55,
    alpha = 0.9,
    outlier.shape = NA,
    linewidth = 0.5
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 19,
    size = 2,
    fill = "black",
    color = "black"
  ) +
  scale_fill_manual(values = fill_cols) +
  labs(
    x = NULL,
    y = "Distance to centroid"
  ) +
  
  annotate(
    "segment",
    x = 1, xend = 2,
    y = y_line, yend = y_line,
    linewidth = 0.5
  ) +
  annotate(
    "segment",
    x = 1, xend = 1,
    y = y_line - (y_max - y_min) * 0.01,
    yend = y_line,
    linewidth = 0.5
  ) +
  annotate(
    "segment",
    x = 2, xend = 2,
    y = y_line - (y_max - y_min) * 0.01,
    yend = y_line,
    linewidth = 0.5
  ) +
  
  # p 值
  annotate(
    "text",
    x = 1.5,
    y = y_text,
    label = p_label,
    size = 5
  ) +
  
  scale_y_continuous(limits = custom_ylim) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 28, hjust = 1),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5)
  )

print(p)

# Output
ggsave(
  filename = "Bacteria_beta_dispersion_custom_boxplot.pdf",
  plot = p,
  width = 6,
  height = 5
)




################################################################################################ PD
library(vegan)
library(picante)
library(ape)
library(ggplot2)

otu_t <- t(otu)

#################### PD
tree <- read.tree("Bacteria_otus.nwk")
tree
length(tree$tip.label)

common_taxa <- intersect(colnames(otu_t), tree$tip.label)
otu_t2 <- otu_t[, common_taxa]
tree2 <- drop.tip(tree, setdiff(tree$tip.label, common_taxa))

length(colnames(otu_t2))
length(tree2$tip.label)

# Calculate PD
pd_res <- pd(otu_t2, tree2, include.root = FALSE)
head(pd_res)

## 6. Merge metadata
pd_df <- data.frame(
  SampleID = rownames(pd_res),
  PD = pd_res$PD,
  SR = pd_res$SR
)

alpha_meta <- merge(metadata, pd_df, by = "SampleID")
## CK, TR
alpha_meta$Group <- factor(alpha_meta$Group, levels = c("CK", "TR"))
## Check
head(alpha_meta)

## Output
write.table(alpha_meta,
            "Bacteria_alpha_diversity.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## P value + BH
vars <- c( "PD")
raw_p <- sapply(vars, function(v) {
  dat <- alpha_meta[, c("Group", v)]
  dat <- dat[complete.cases(dat), ]
  wilcox.test(dat[[v]] ~ dat$Group, exact = FALSE)$p.value
})
p_adj <- p.adjust(raw_p, method = "BH")
stat_res <- data.frame(
  Index = vars,
  P_raw = raw_p,
  P_adj = p_adj,
  stringsAsFactors = FALSE
)
print(stat_res)

## Output
write.table(stat_res, "Bacteria_alpha_Wilcoxon_BH.txt", sep = "\t", quote = FALSE, row.names = FALSE)

##  Color
my_fill <- c("CK" = "#EFB44F", "TR" = "#5E97C9")

plot_alpha <- function(data, yvar, stat_table, ylab = NULL,
                       out_pdf = NULL, out_png = NULL) {
  p_use <- stat_table$P_adj[stat_table$Index == yvar]
  p_lab <- paste0("Wilcoxon, BH-adjusted p = ", signif(p_use, 3))
  dat <- data[, c("Group", yvar)]
  dat <- dat[complete.cases(dat), ]
  colnames(dat)[2] <- "Value"
  y_max <- max(dat$Value, na.rm = TRUE)
  y_min <- min(dat$Value, na.rm = TRUE)
  y_pos <- y_max + (y_max - y_min) * 0.08
  p <- ggplot(dat, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot(width = 0.55, alpha = 0.9, outlier.shape = NA, linewidth = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2,fill = "black", color = "black") +
    scale_fill_manual(values = my_fill) +
    labs(x = NULL, y = ifelse(is.null(ylab), yvar, ylab)) +
    annotate("text", x = 1.5, y = y_pos, label = p_lab, size = 5) +
    theme_classic(base_size = 16) +
    theme(  legend.position = "none",  axis.line = element_line(linewidth = 0.5), axis.ticks = element_line(linewidth = 0.5),
            axis.text.x = element_text(angle = 25, hjust = 1), plot.margin = margin(10, 15, 10, 10)) +
    coord_cartesian(ylim = c(0, 350))
  print(p)
  if (!is.null(out_pdf)) {
    ggsave(out_pdf, p, width = 6, height = 5)
  }
  return(p)
}

p3 <- plot_alpha(alpha_meta,
                 yvar = "PD",
                 stat_table = stat_res,
                 ylab = "PD",
                 out_pdf = "Bacteria_PD_boxplot.pdf")



##################################################################################################### Niche Overlap

library(parallel)
library(dplyr)
library(ggplot2)
library(openxlsx)

start_time <- Sys.time()

otu_mat <- as.matrix(otu)
mode(otu_mat) <- "numeric"
asv_rel <- sweep(otu_mat, 2, colSums(otu_mat), "/")

sample_names <- colnames(asv_rel)
pairs <- t(combn(sample_names, 2))

schoener <- function(p, q) {
  1 - 0.5 * sum(abs(p - q))
}

cl <- makeCluster(24)
clusterExport(cl, varlist = c("asv_rel", "schoener"))

schoener_values <- parApply(cl, pairs, 1, function(x) {
  s1 <- x[1]
  s2 <- x[2]
  schoener(asv_rel[, s1], asv_rel[, s2])
})

stopCluster(cl)

schoener_df <- data.frame(
  Sample1 = pairs[, 1],
  Sample2 = pairs[, 2],
  Schoener = schoener_values,
  stringsAsFactors = FALSE
)

metadata$SampleID <- trimws(metadata$SampleID)
metadata$Group <- trimws(metadata$Group)

metadata2 <- metadata %>%
  filter(SampleID %in% sample_names) %>%
  select(SampleID, Group)

schoener_df <- schoener_df %>%
  left_join(metadata2, by = c("Sample1" = "SampleID")) %>%
  rename(Group1 = Group) %>%
  left_join(metadata2, by = c("Sample2" = "SampleID")) %>%
  rename(Group2 = Group)

group_within <- schoener_df %>%
  filter(Group1 == Group2) %>%
  mutate(Group = factor(Group1, levels = c("CK", "TR")))

summary_df <- group_within %>%
  group_by(Group) %>%
  summarise( mean = mean(Schoener, na.rm = TRUE), sd   = sd(Schoener, na.rm = TRUE), n    = n(), se   = sd / sqrt(n),  .groups = "drop" )
print(summary_df)

wilcox_res <- wilcox.test(Schoener ~ Group, data = group_within, exact = FALSE, p.adjust.method = "BH")
print(wilcox_res)
p_value <- wilcox_res$p.value
p_label <- paste0("Wilcoxon p = ", signif(p_value, 3))

write.xlsx(summary_df, "Bacteria_Schoener_Index_Summary.xlsx", rowNames = FALSE)
write.xlsx(schoener_df, "Bacteria_Schoener_Index_All_Pairs.xlsx", rowNames = FALSE)
write.xlsx(group_within, "Bacteria_Schoener_Index_Within_Group_Pairs.xlsx", rowNames = FALSE)

fill_cols <- c("CK" = "#F7AF34", "TR" = "#448DCD")

y_min <- min(group_within$Schoener, na.rm = TRUE)
y_max <- max(group_within$Schoener, na.rm = TRUE)

y_line <- y_max * 0.23
y_text <- y_max * 0.25

p <- ggplot(group_within, aes(x = Group, y = Schoener, fill = Group)) +
  geom_boxplot(
    width = 0.55, alpha = 0.9, outlier.shape = NA, linewidth = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2,fill = "black", color = "black") +
  scale_fill_manual(values = fill_cols) +
  labs(x = NULL, y = "Schoener's index (within-group)") +

  annotate("segment", x = 1, xend = 2, y = y_line, yend = y_line, linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, y = y_line * 0.995, yend = y_line, linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, y = y_line * 0.995, yend = y_line, linewidth = 0.5) +

  annotate("text", x = 1.5, y = y_text, label = p_label, size = 5) +
  scale_y_continuous(limits = c(y_min, y_text * 1.03)) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 28, hjust = 1),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5)
  )
print(p)
### output
ggsave("Bacteria_Schoener_within_groups_boxplot.pdf", p, width = 6, height = 5)

end_time <- Sys.time()
print(end_time - start_time)






##################################################################################################### ALDEx2
library(ALDEx2)
library(dplyr)
library(ggplot2)
library(ggrepel)

asv_table <- read.delim("Bacteria_otutab_rare.txt",
                        row.names = 1,
                        check.names = FALSE,
                        stringsAsFactors = FALSE)

taxonomy <- read.delim("Bacteria_taxonomy.txt",
                       row.names = 1,
                       check.names = FALSE,
                       stringsAsFactors = FALSE)

metadata <- read.delim("Bacteria_metadata_raw.txt",
                       check.names = FALSE,
                       stringsAsFactors = FALSE)

common_samples <- intersect(metadata$SampleID, colnames(asv_table))

asv_table <- asv_table[, common_samples, drop = FALSE]
metadata <- metadata[match(colnames(asv_table), metadata$SampleID), ]

stopifnot(all(metadata$SampleID == colnames(asv_table)))

asv_table[] <- lapply(asv_table, function(x) as.numeric(as.character(x)))
asv_table[is.na(asv_table)] <- 0

asv_table$OTUID <- rownames(asv_table)
taxonomy$OTUID <- rownames(taxonomy)
asv_taxa <- left_join(asv_table, taxonomy, by = "OTUID")

run_aldex_on_level <- function(taxa_level) {
  message("Running ALDEx2 on level: ", taxa_level)
  remove_cols <- c("OTUID", "Kingdom", "Phylum", "Class",
                   "Order", "Family", "Genus", "Species")
  level_vec <- asv_taxa[[taxa_level]]
  level_vec[is.na(level_vec) | level_vec == ""] <- paste0("Unclassified_", taxa_level)
  taxa_data <- asv_taxa %>%
    dplyr::select(-dplyr::any_of(remove_cols)) %>%
    bind_cols(Level = level_vec) %>%
    group_by(Level) %>%
    summarise(across(everything(), ~ sum(as.numeric(.x), na.rm = TRUE)),
              .groups = "drop")
  taxa_df <- as.data.frame(taxa_data)
  rownames(taxa_df) <- taxa_df$Level
  taxa_df$Level <- NULL
  taxa_df <- taxa_df[, metadata$SampleID, drop = FALSE]
  taxa_df <- taxa_df[rowSums(taxa_df) > 0, , drop = FALSE]
  taxa_mat <- as.matrix(taxa_df)
  storage.mode(taxa_mat) <- "numeric"
  group <- as.character(metadata$Group)
  set.seed(1)
  aldex_clr <- aldex.clr(taxa_mat,
                         conds = group,
                         mc.samples = 128,
                         denom = "all",
                         verbose = TRUE)
  aldex_res <- aldex.ttest(aldex_clr, paired.test = FALSE)
  aldex_eff <- aldex.effect(aldex_clr)
  res <- cbind(aldex_res, aldex_eff)
  res$Taxon <- rownames(res)
  res <- res[, c("Taxon",
                 "we.ep", "we.eBH",
                 "wi.ep", "wi.eBH",
                 "effect",
                 "diff.btw",
                 "diff.win")]
  res <- res[order(res$wi.eBH), ]
  res$Significance <- ifelse(res$wi.eBH < 0.001, "***",
                             ifelse(res$wi.eBH < 0.01, "**",
                                    ifelse(res$wi.eBH < 0.05, "*", "")))
  write.csv(res,
            paste0("Bacteria_", taxa_level, "_aldex_result.csv"),
            row.names = FALSE)
  return(list(res = res, taxa_df = taxa_df))
}


############################# plot, top10
plot_top10_effect <- function(res, taxa_df, level_name) {
  
  library(ggplot2)
  
  ## 1. Remove Unassigned / Unclassified
  remove_pattern <- "unassigned|uncultured|unclassified|unknown|NA"
  
  keep_taxa <- !grepl(remove_pattern,
                      rownames(taxa_df),
                      ignore.case = TRUE)
  
  taxa_df_plot <- taxa_df[keep_taxa, , drop = FALSE]
  
  ## 2. Top10
  rel_abund <- sweep(taxa_df_plot, 2, colSums(taxa_df_plot), "/")
  rel_abund[is.na(rel_abund)] <- 0
  
  mean_abund <- rowMeans(rel_abund)
  
  top10_taxa <- names(sort(mean_abund, decreasing = TRUE))[1:min(10, length(mean_abund))]
  
  ## 3.  Top10_ALDEx2
  res_plot <- res[res$Taxon %in% top10_taxa, , drop = FALSE]
  
  res_plot <- res_plot[order(res_plot$effect), ]
  res_plot$Taxon <- factor(res_plot$Taxon, levels = res_plot$Taxon)

  ## 4. Plot
  p <- ggplot(res_plot, aes(x = effect, y = Taxon)) +
    geom_col(aes(fill = effect > 0), width = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.6) +
    
    geom_text(aes(label = Significance),
              hjust = ifelse(res_plot$effect > 0, -0.3, 1.3),
              size = 5) +
    scale_fill_manual(values = c("TRUE" = "#F7AF34",
                                 "FALSE" = "#448DCD"),
                      guide = "none") +
    
    labs(x = "Effect size of difference",
         y = NULL,
         title = paste0("Top 10 ", level_name)) +
    
    theme_classic(base_size = 14) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(10, 35, 10, 10))
  ## 5. Output
  ggsave(paste0("Bacteria_", level_name, "_Top10_effect_noUnassigned.pdf"),
         p,
         width = 7,
         height = 5,
         device = cairo_pdf)
  
  print(p)
}

## Phylum
phylum_out <- run_aldex_on_level("Phylum")
plot_top10_effect(phylum_out$res,
                  phylum_out$taxa_df,
                  "Phylum")
## Class
class_out <- run_aldex_on_level("Class")
plot_top10_effect(class_out$res,
                  class_out$taxa_df,
                  "Class")



################################################################################## Separate ASV for network analysis
library(dplyr)

metadata$SampleID   <- trimws(metadata$SampleID)
metadata$Group      <- trimws(metadata$Group)
colnames(otu)   <- trimws(colnames(otu))
rownames(otu)   <- trimws(rownames(otu))


overlap_cols <- length(intersect(metadata$SampleID, colnames(otu)))
overlap_rows <- length(intersect(metadata$SampleID, rownames(otu)))

if (is.na(overlap_cols)) overlap_cols <- 0
if (is.na(overlap_rows)) overlap_rows <- 0

if (overlap_rows > overlap_cols) {
  message("检测到样本名在行名中，已自动转置 OTU 表（现在列=样本）。")
  otu <- t(otu)
}

common_samples <- intersect(metadata$SampleID, colnames(otu))
missing_in_otu  <- setdiff(metadata$SampleID, colnames(otu))
missing_in_meta <- setdiff(colnames(otu), metadata$SampleID)

message(sprintf("共有样本匹配: %d", length(common_samples)))
if (length(missing_in_otu) > 0)
  message(sprintf("metadata 有但 OTU 表没有的样本（前10个）：%s",
                  paste(head(missing_in_otu, 10), collapse = ", ")))
if (length(missing_in_meta) > 0)
  message(sprintf("OTU 表有但 metadata 没有的样本（前10个）：%s",
                  paste(head(missing_in_meta, 10), collapse = ", ")))

meta_use <- metadata %>% filter(SampleID %in% common_samples)
otu_use  <- otu[, meta_use$SampleID, drop = FALSE]

otu_split <- split(meta_use$SampleID, meta_use$Group)
otu_split <- lapply(otu_split, function(smp) otu_use[, smp, drop = FALSE])

for (g in names(otu_split)) {
  fn <- paste0("Bacteria_otu_group_", g, ".txt")
  write.table(otu_split[[g]], file = fn, sep = "\t", quote = FALSE, col.names = NA)
}


###############################################################################################  Control Network
suppressPackageStartupMessages({
  library(Hmisc)
  library(igraph)
})

.trim <- function(x) gsub("\\s+$","", gsub("^\\s+","", as.character(x)))
.safe_read_table <- function(path){
  if (grepl("\\.csv$", path, ignore.case=TRUE)) {
    read.csv(path, header=TRUE, row.names=1, check.names=FALSE)
  } else {
    read.table(path, header=TRUE, row.names=1, sep="\t",
               check.names=FALSE, fill=TRUE, quote="", comment.char="")
  }
}

build_network <- function(
    otu_input,                 
    tax_input = NULL,        
    pre_filter = TRUE,         
    prev_min = 0.20,         
    mean_rel_min = 1e-4,       
    pool_global = FALSE,       
    prevalence_prop = 0,  
    total_abund_prop = 0,  
    rho_thr = 0.6,
    q_thr = 0.01,
    allow_negative = TRUE,
    # Output
    output_dir = "Bacteria_CK_network_out",
    prefix = "Bacteria_CK"
){
  # Input OTU
  otu0 <- if (is.character(otu_input)) .safe_read_table(otu_input) else as.data.frame(otu_input, check.names=FALSE)
  stopifnot(!is.null(rownames(otu0)))
  rn <- rownames(otu0); cn <- colnames(otu0)
  asv_in_rows <- mean(grepl("^(ASV|OTU)_|[A-Za-z]", rn)) >=
    mean(grepl("^(ASV|OTU)_|[A-Za-z]", cn))
  otu <- if (asv_in_rows) t(otu0) else otu0
  otu <- otu[rowSums(otu) > 0, , drop=FALSE]
  otu <- otu[, colSums(otu) > 0, drop=FALSE]
  if (nrow(otu) < 2 || ncol(otu) < 2) stop("有效行/列不足：请检查 OTU 表。")
  rs <- rowSums(otu); rs[rs==0] <- 1
  otu_rel <- sweep(otu, 1, rs, "/")

  if (isTRUE(pre_filter)) {
    prev <- colMeans(otu_rel > 0)
    if (isTRUE(pool_global)) {
      abund <- colSums(otu) / sum(colSums(otu))
    } else {
      abund <- colMeans(otu_rel) 
    }
    keep <- (prev >= prev_min) & (abund >= mean_rel_min)
    otu_rel <- otu_rel[, keep, drop=FALSE]
    if (ncol(otu_rel) < 2) stop("预过滤后 ASV < 2，请调低 prev_min / mean_rel_min。")
  }
  
  if (prevalence_prop > 0 || total_abund_prop > 0) {
    prev_prop   <- colMeans(otu_rel > 0)
    total_abund <- colMeans(otu_rel)  
    keep2 <- (prev_prop >= prevalence_prop) & (total_abund >= total_abund_prop)
    otu_rel <- otu_rel[, keep2, drop=FALSE]
    if (ncol(otu_rel) < 2) stop("叠加过滤后 ASV < 2，请调低 prevalence_prop / total_abund_prop。")
  }
  
  #  Spearman  + FDR 
  rc <- Hmisc::rcorr(as.matrix(otu_rel), type="spearman")
  R <- rc$r; P <- rc$P
  ut <- upper.tri(R)
  cor_df <- data.frame(
    source   = rownames(R)[col(R)[ut]],
    target   = rownames(R)[row(R)[ut]],
    rho      = R[ut],
    weight   = abs(R[ut]),
    p        = P[ut],
    q        = p.adjust(P[ut], method="BH"),
    cor_type = ifelse(R[ut] > 0, "Positive", "Negative"),
    stringsAsFactors = FALSE
  )
  sig <- if (allow_negative) subset(cor_df, weight >= rho_thr & q < q_thr)
  else                 subset(cor_df, rho    >= rho_thr & q < q_thr)
  if (nrow(sig) == 0) stop("筛选后无边；请放宽 rho_thr 或 q_thr。")
  
  # Construct network
  g <- igraph::graph_from_data_frame(sig, directed=FALSE)
  g <- delete_vertices(g, V(g)[degree(g) == 0])
  
  # taxonomy
  if (!is.null(tax_input)) {
    tax0 <- if (is.character(tax_input)) {
      if (grepl("\\.csv$", tax_input, ignore.case=TRUE)) {
        read.csv(tax_input, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
      } else {
        read.table(tax_input, header=TRUE, sep="\t", check.names=FALSE,
                   stringsAsFactors=FALSE, quote="", comment.char="", fill=TRUE)
      }
    } else as.data.frame(tax_input, check.names=FALSE)
    V(g)$name <- .trim(V(g)$name)
    if (!is.null(rownames(tax0)) && !all(rownames(tax0) == as.character(seq_len(nrow(tax0))))) {
      tax0$.__rowid__ <- .trim(rownames(tax0))
    }
    id_candidates <- c("OTUID","ASV","FeatureID","ID","#OTU ID")
    idcol <- id_candidates[id_candidates %in% colnames(tax0)]
    if (length(idcol) == 0) idcol <- if (!is.null(tax0$.__rowid__)) ".__rowid__" else colnames(tax0)[1]
    tax0[[idcol[1]]] <- .trim(tax0[[idcol[1]]])
    tax0 <- tax0[!duplicated(tax0[[idcol[1]]]), , drop=FALSE]
    rownames(tax0) <- tax0[[idcol[1]]]
    idx <- match(V(g)$name, rownames(tax0))
    message("taxonomy 匹配到的节点数: ", sum(!is.na(idx)), " / ", vcount(g))
    tax_cols <- intersect(colnames(tax0), c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
    if (length(tax_cols) > 0 && any(!is.na(idx))) {
      for (nm in tax_cols) {
        vals <- rep(NA_character_, vcount(g)); ok <- !is.na(idx)
        vals[ok] <- as.character(tax0[[nm]][idx[ok]])
        g <- igraph::set_vertex_attr(g, nm, value = vals)
      }
    } else {
      warning("taxonomy 未匹配到常见分类列或匹配数为 0，跳过写入。")
    }
  }
  
  # Module + Zi/Pi + Keystone
  cl <- cluster_louvain(g)
  V(g)$module <- as.character(membership(cl))
  deg_total <- degree(g)
  modules <- V(g)$module; nodes <- V(g)$name
  mods_unique <- unique(modules)
  adj_list <- adjacent_vertices(g, V(g))
  
  k_is_self <- numeric(vcount(g))
  k_is_by_mod <- vector("list", length=vcount(g))
  for (i in seq_along(nodes)) {
    nbrs <- names(adj_list[[i]])
    nbr_mods <- modules[match(nbrs, nodes)]
    tab <- table(nbr_mods)
    ki_vec <- setNames(numeric(length(mods_unique)), mods_unique)
    ki_vec[names(tab)] <- as.numeric(tab)
    k_is_by_mod[[i]] <- ki_vec
    k_is_self[i] <- ki_vec[modules[i]]
  }
  kstats <- lapply(mods_unique, function(m){
    idxm <- which(modules == m); kis <- k_is_self[idxm]
    c(mean = mean(kis), sd = stats::sd(kis))
  }); names(kstats) <- mods_unique
  
  Zi <- sapply(seq_along(nodes), function(i){
    m <- modules[i]; mu <- kstats[[m]]["mean"]; sdv <- kstats[[m]]["sd"]
    ifelse(is.na(sdv) || sdv == 0, 0, (k_is_self[i] - mu)/sdv)
  })
  Pi <- sapply(seq_along(nodes), function(i){
    ki <- deg_total[i]; if (ki == 0) 0 else 1 - sum((k_is_by_mod[[i]]/ki)^2, na.rm=TRUE)
  })
  V(g)$Zi <- Zi; V(g)$Pi <- Pi
  V(g)$role <- ifelse(Zi >= 2.5 & Pi < 0.62, "module_hub",
                      ifelse(Zi < 2.5 & Pi >= 0.62, "connector",
                             ifelse(Zi >= 2.5 & Pi >= 0.62, "network_hub","peripheral")))
  
  # Output
  pos_num <- sum(E(g)$rho > 0); neg_num <- sum(E(g)$rho < 0)
  stats <- data.frame(
    nodes_num = vcount(g), edges_num = ecount(g),
    positive.corr_num  = pos_num, positive.corr_prop = pos_num/ecount(g),
    negative.corr_num  = neg_num, negative.corr_prop = neg_num/ecount(g),
    average_degree = mean(degree(g)),
    average_shortest_path_length = suppressWarnings(mean_distance(g, directed=FALSE)),
    network_diameter = suppressWarnings(diameter(g, directed=FALSE)),
    network_density  = edge_density(g),
    clustering_coefficient = transitivity(g, type="global"),
    modularity = modularity(g, membership(cl), directed=FALSE)
  )
  
  dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  f_nodes <- file.path(output_dir, sprintf("%s_nodes_%s.csv",  prefix, ts))
  f_edges <- file.path(output_dir, sprintf("%s_edges_%s.csv",  prefix, ts))
  f_stats <- file.path(output_dir, sprintf("%s_stats_%s.csv",  prefix, ts))
  f_graph <- file.path(output_dir, sprintf("%s_graph_%s.graphml", prefix, ts))
  
  nodes_df <- data.frame(
    id=V(g)$name, module=V(g)$module, degree=deg_total,
    Zi=V(g)$Zi, Pi=V(g)$Pi, role=V(g)$role, stringsAsFactors=FALSE
  )
  keep_attr <- setdiff(vertex_attr_names(g), c("name","module","Zi","Pi","role"))
  for (att in keep_attr) nodes_df[[att]] <- igraph::get.vertex.attribute(g, att)
  
  edges_df <- igraph::as_data_frame(g, what="edges")
  colnames(edges_df)[1:2] <- c("source","target")
  
  write.csv(nodes_df, f_nodes, row.names=FALSE)
  write.csv(edges_df, f_edges, row.names=FALSE)
  write.csv(stats,    f_stats, row.names=FALSE)
  igraph::write_graph(g, f_graph, format="graphml")
  
  list(graph=g, nodes=nodes_df, edges=edges_df, stats=stats,
       files=list(nodes=f_nodes, edges=f_edges, stats=f_stats, graphml=f_graph))
}

otu_file <- "Bacteria_otu_group_CK.txt"
tax_file <- "Bacteria_taxonomy.txt"

otu_raw <- read.table(otu_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
otu_input_ok <- t(otu_raw)

# Construct
res <- build_network(
  otu_input = otu_file,
  tax_input = tax_file,
  pre_filter = TRUE,  
  prev_min = 0.20,
  mean_rel_min = 1e-4,
  pool_global = FALSE, 
  prevalence_prop = 0, total_abund_prop = 0,  
  rho_thr = 0.6, q_thr = 0.01, allow_negative = TRUE,
  output_dir = "Bacteria_CK_network_out", prefix = "Bacteria_CK"
)
# Check results
res$stats
head(res$edges) 
head(res$nodes)   
#   nodes_num edges_num positive.corr_num positive.corr_prop negative.corr_num negative.corr_prop average_degree
# 1       507      5519              5516          0.9994564                 3       0.0005435767        21.7712
#   average_shortest_path_length network_diameter network_density clustering_coefficient modularity
# 1                     2.867545         8.187147      0.04302609              0.7047182  0.5539679




###############################################################################################  Rotation Network
suppressPackageStartupMessages({
  library(Hmisc)
  library(igraph)
})

.trim <- function(x) gsub("\\s+$","", gsub("^\\s+","", as.character(x)))
.safe_read_table <- function(path){
  if (grepl("\\.csv$", path, ignore.case=TRUE)) {
    read.csv(path, header=TRUE, row.names=1, check.names=FALSE)
  } else {
    read.table(path, header=TRUE, row.names=1, sep="\t",
               check.names=FALSE, fill=TRUE, quote="", comment.char="")
  }
}

build_network <- function(
    otu_input,                 
    tax_input = NULL,          
    pre_filter = TRUE,         
    prev_min = 0.20,          
    mean_rel_min = 1e-4,       
    pool_global = FALSE,       
    prevalence_prop = 0,       
    total_abund_prop = 0,   
    rho_thr = 0.6,
    q_thr = 0.01,
    allow_negative = TRUE,
    output_dir = "Bacteria_TR_network_out",
    prefix = "Bacteria_TR"
){

  otu0 <- if (is.character(otu_input)) .safe_read_table(otu_input) else as.data.frame(otu_input, check.names=FALSE)
  stopifnot(!is.null(rownames(otu0)))
  rn <- rownames(otu0); cn <- colnames(otu0)

  asv_in_rows <- mean(grepl("^(ASV|OTU)_|[A-Za-z]", rn)) >=
    mean(grepl("^(ASV|OTU)_|[A-Za-z]", cn))
  otu <- if (asv_in_rows) t(otu0) else otu0
  otu <- otu[rowSums(otu) > 0, , drop=FALSE]
  otu <- otu[, colSums(otu) > 0, drop=FALSE]
  if (nrow(otu) < 2 || ncol(otu) < 2) stop("有效行/列不足：请检查 OTU 表。")

  rs <- rowSums(otu); rs[rs==0] <- 1
  otu_rel <- sweep(otu, 1, rs, "/")

  if (isTRUE(pre_filter)) {
    prev <- colMeans(otu_rel > 0)
    if (isTRUE(pool_global)) {
      abund <- colSums(otu) / sum(colSums(otu))
    } else {
      abund <- colMeans(otu_rel)
    }
    keep <- (prev >= prev_min) & (abund >= mean_rel_min)
    otu_rel <- otu_rel[, keep, drop=FALSE]
    if (ncol(otu_rel) < 2) stop("预过滤后 ASV < 2，请调低 prev_min / mean_rel_min。")
  }
  

  if (prevalence_prop > 0 || total_abund_prop > 0) {
    prev_prop   <- colMeans(otu_rel > 0)
    total_abund <- colMeans(otu_rel)  # 已是相对丰度矩阵，用均值代表总体强度
    keep2 <- (prev_prop >= prevalence_prop) & (total_abund >= total_abund_prop)
    otu_rel <- otu_rel[, keep2, drop=FALSE]
    if (ncol(otu_rel) < 2) stop("叠加过滤后 ASV < 2，请调低 prevalence_prop / total_abund_prop。")
  }
  
  #  Spearman + FDR 
  rc <- Hmisc::rcorr(as.matrix(otu_rel), type="spearman")
  R <- rc$r; P <- rc$P
  ut <- upper.tri(R)
  cor_df <- data.frame(
    source   = rownames(R)[col(R)[ut]],
    target   = rownames(R)[row(R)[ut]],
    rho      = R[ut],
    weight   = abs(R[ut]),
    p        = P[ut],
    q        = p.adjust(P[ut], method="BH"),
    cor_type = ifelse(R[ut] > 0, "Positive", "Negative"),
    stringsAsFactors = FALSE
  )
  sig <- if (allow_negative) subset(cor_df, weight >= rho_thr & q < q_thr)
  else                 subset(cor_df, rho    >= rho_thr & q < q_thr)
  if (nrow(sig) == 0) stop("筛选后无边；请放宽 rho_thr 或 q_thr。")
  
  g <- igraph::graph_from_data_frame(sig, directed=FALSE)
  g <- delete_vertices(g, V(g)[degree(g) == 0])
  
  if (!is.null(tax_input)) {
    tax0 <- if (is.character(tax_input)) {
      if (grepl("\\.csv$", tax_input, ignore.case=TRUE)) {
        read.csv(tax_input, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
      } else {
        read.table(tax_input, header=TRUE, sep="\t", check.names=FALSE,
                   stringsAsFactors=FALSE, quote="", comment.char="", fill=TRUE)
      }
    } else as.data.frame(tax_input, check.names=FALSE)
    V(g)$name <- .trim(V(g)$name)
    if (!is.null(rownames(tax0)) && !all(rownames(tax0) == as.character(seq_len(nrow(tax0))))) {
      tax0$.__rowid__ <- .trim(rownames(tax0))
    }
    id_candidates <- c("OTUID","ASV","FeatureID","ID","#OTU ID")
    idcol <- id_candidates[id_candidates %in% colnames(tax0)]
    if (length(idcol) == 0) idcol <- if (!is.null(tax0$.__rowid__)) ".__rowid__" else colnames(tax0)[1]
    tax0[[idcol[1]]] <- .trim(tax0[[idcol[1]]])
    tax0 <- tax0[!duplicated(tax0[[idcol[1]]]), , drop=FALSE]
    rownames(tax0) <- tax0[[idcol[1]]]
    idx <- match(V(g)$name, rownames(tax0))
    message("taxonomy 匹配到的节点数: ", sum(!is.na(idx)), " / ", vcount(g))
    tax_cols <- intersect(colnames(tax0), c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
    if (length(tax_cols) > 0 && any(!is.na(idx))) {
      for (nm in tax_cols) {
        vals <- rep(NA_character_, vcount(g)); ok <- !is.na(idx)
        vals[ok] <- as.character(tax0[[nm]][idx[ok]])
        g <- igraph::set_vertex_attr(g, nm, value = vals)
      }
    } else {
      warning("taxonomy 未匹配到常见分类列或匹配数为 0，跳过写入。")
    }
  }
  
  cl <- cluster_louvain(g)
  V(g)$module <- as.character(membership(cl))
  deg_total <- degree(g)
  modules <- V(g)$module; nodes <- V(g)$name
  mods_unique <- unique(modules)
  adj_list <- adjacent_vertices(g, V(g))
  
  k_is_self <- numeric(vcount(g))
  k_is_by_mod <- vector("list", length=vcount(g))
  for (i in seq_along(nodes)) {
    nbrs <- names(adj_list[[i]])
    nbr_mods <- modules[match(nbrs, nodes)]
    tab <- table(nbr_mods)
    ki_vec <- setNames(numeric(length(mods_unique)), mods_unique)
    ki_vec[names(tab)] <- as.numeric(tab)
    k_is_by_mod[[i]] <- ki_vec
    k_is_self[i] <- ki_vec[modules[i]]
  }
  kstats <- lapply(mods_unique, function(m){
    idxm <- which(modules == m); kis <- k_is_self[idxm]
    c(mean = mean(kis), sd = stats::sd(kis))
  }); names(kstats) <- mods_unique
  
  Zi <- sapply(seq_along(nodes), function(i){
    m <- modules[i]; mu <- kstats[[m]]["mean"]; sdv <- kstats[[m]]["sd"]
    ifelse(is.na(sdv) || sdv == 0, 0, (k_is_self[i] - mu)/sdv)
  })
  Pi <- sapply(seq_along(nodes), function(i){
    ki <- deg_total[i]; if (ki == 0) 0 else 1 - sum((k_is_by_mod[[i]]/ki)^2, na.rm=TRUE)
  })
  V(g)$Zi <- Zi; V(g)$Pi <- Pi
  V(g)$role <- ifelse(Zi >= 2.5 & Pi < 0.62, "module_hub",
                      ifelse(Zi < 2.5 & Pi >= 0.62, "connector",
                             ifelse(Zi >= 2.5 & Pi >= 0.62, "network_hub","peripheral")))

  pos_num <- sum(E(g)$rho > 0); neg_num <- sum(E(g)$rho < 0)
  stats <- data.frame(
    nodes_num = vcount(g), edges_num = ecount(g),
    positive.corr_num  = pos_num, positive.corr_prop = pos_num/ecount(g),
    negative.corr_num  = neg_num, negative.corr_prop = neg_num/ecount(g),
    average_degree = mean(degree(g)),
    average_shortest_path_length = suppressWarnings(mean_distance(g, directed=FALSE)),
    network_diameter = suppressWarnings(diameter(g, directed=FALSE)),
    network_density  = edge_density(g),
    clustering_coefficient = transitivity(g, type="global"),
    modularity = modularity(g, membership(cl), directed=FALSE)
  )
  
  dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  f_nodes <- file.path(output_dir, sprintf("%s_nodes_%s.csv",  prefix, ts))
  f_edges <- file.path(output_dir, sprintf("%s_edges_%s.csv",  prefix, ts))
  f_stats <- file.path(output_dir, sprintf("%s_stats_%s.csv",  prefix, ts))
  f_graph <- file.path(output_dir, sprintf("%s_graph_%s.graphml", prefix, ts))
  
  nodes_df <- data.frame(
    id=V(g)$name, module=V(g)$module, degree=deg_total,
    Zi=V(g)$Zi, Pi=V(g)$Pi, role=V(g)$role, stringsAsFactors=FALSE
  )
  keep_attr <- setdiff(vertex_attr_names(g), c("name","module","Zi","Pi","role"))
  for (att in keep_attr) nodes_df[[att]] <- igraph::get.vertex.attribute(g, att)
  
  edges_df <- igraph::as_data_frame(g, what="edges")
  colnames(edges_df)[1:2] <- c("source","target")
  
  write.csv(nodes_df, f_nodes, row.names=FALSE)
  write.csv(edges_df, f_edges, row.names=FALSE)
  write.csv(stats,    f_stats, row.names=FALSE)
  igraph::write_graph(g, f_graph, format="graphml")
  
  list(graph=g, nodes=nodes_df, edges=edges_df, stats=stats,
       files=list(nodes=f_nodes, edges=f_edges, stats=f_stats, graphml=f_graph))
}

otu_file <- "Bacteria_otu_group_TR.txt"
tax_file <- "Bacteria_taxonomy.txt"

otu_raw <- read.table(otu_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
otu_input_ok <- t(otu_raw)

res <- build_network(
  otu_input = otu_file,
  tax_input = tax_file,
  pre_filter = TRUE,  
  prev_min = 0.20,
  mean_rel_min = 1e-4,
  pool_global = FALSE, 
  prevalence_prop = 0, total_abund_prop = 0,  
  rho_thr = 0.6, q_thr = 0.01, allow_negative = TRUE,
  output_dir = "Bacteria_TR_network_out", prefix = "Bacteria_TR"
)

res$stats
head(res$edges)    
head(res$nodes)   
#   nodes_num edges_num positive.corr_num positive.corr_prop negative.corr_num negative.corr_prop average_degree
# 1       506     13445             13444          0.9999256                 1       7.437709e-05       53.14229
#   average_shortest_path_length network_diameter network_density clustering_coefficient modularity
# 1                     2.687809         11.19823       0.1052323              0.7115469  0.1895581





########################################################################################################## Subnetwork

# Full pipeline:
# 1) Build microbial co-occurrence network (Spearman + BH-FDR)
# 2) Extract sample-specific subnetworks from the global network
# 3) Compute topology metrics per sample for soil correlation
#    (sample metrics step is parallelized with workers)

suppressPackageStartupMessages({
  library(Hmisc)
  library(igraph)
  library(dplyr)
  library(future)
  library(future.apply)
})

# Utils
.trim <- function(x) gsub("\\s+$","", gsub("^\\s+","", as.character(x)))

.safe_read_table <- function(path){
  if (grepl("\\.csv$", path, ignore.case=TRUE)) {
    read.csv(path, header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
  } else {
    read.table(path, header=TRUE, row.names=1, sep="\t",
               check.names=FALSE, fill=TRUE, quote="", comment.char="",
               stringsAsFactors=FALSE)
  }
}

# Build global network (Spearman + BH-FDR)
build_network <- function(
    otu_input,            
    tax_input = NULL,          
    pre_filter = TRUE,
    prev_min = 0.20,
    mean_rel_min = 1e-4,
    pool_global = FALSE,
    prevalence_prop = 0,
    total_abund_prop = 0,
    rho_thr = 0.6,
    q_thr = 0.01,
    allow_negative = TRUE,
    # 导出
    output_dir = "Bacteria_Subnetwork_out",
    prefix = "Bacteria_Subnetwork"
){
  otu0 <- if (is.character(otu_input)) .safe_read_table(otu_input) else as.data.frame(otu_input, check.names=FALSE)
  stopifnot(!is.null(rownames(otu0)))
  
  rn <- rownames(otu0); cn <- colnames(otu0)

  asv_in_rows <- mean(grepl("^(ASV|OTU)_|[A-Za-z]", rn)) >=
    mean(grepl("^(ASV|OTU)_|[A-Za-z]", cn))
  otu <- if (asv_in_rows) t(otu0) else otu0
  
  otu <- otu[rowSums(otu) > 0, , drop=FALSE]
  otu <- otu[, colSums(otu) > 0, drop=FALSE]
  if (nrow(otu) < 2 || ncol(otu) < 2) stop("有效行/列不足：请检查 OTU 表。")
  
  rs <- rowSums(otu); rs[rs==0] <- 1
  otu_rel <- sweep(otu, 1, rs, "/")

  if (isTRUE(pre_filter)) {
    prev <- colMeans(otu_rel > 0)
    abund <- if (isTRUE(pool_global)) colSums(otu) / sum(colSums(otu)) else colMeans(otu_rel)
    keep <- (prev >= prev_min) & (abund >= mean_rel_min)
    otu_rel <- otu_rel[, keep, drop=FALSE]
    if (ncol(otu_rel) < 2) stop("预过滤后 ASV < 2，请调低 prev_min / mean_rel_min。")
  }

  if (prevalence_prop > 0 || total_abund_prop > 0) {
    prev_prop   <- colMeans(otu_rel > 0)
    total_abund <- colMeans(otu_rel)
    keep2 <- (prev_prop >= prevalence_prop) & (total_abund >= total_abund_prop)
    otu_rel <- otu_rel[, keep2, drop=FALSE]
    if (ncol(otu_rel) < 2) stop("叠加过滤后 ASV < 2，请调低 prevalence_prop / total_abund_prop。")
  }

  rc <- Hmisc::rcorr(as.matrix(otu_rel), type="spearman")
  R <- rc$r; P <- rc$P
  ut <- upper.tri(R)
  
  cor_df <- data.frame(
    source   = rownames(R)[col(R)[ut]],
    target   = rownames(R)[row(R)[ut]],
    rho      = R[ut],
    weight   = abs(R[ut]),
    p        = P[ut],
    q        = p.adjust(P[ut], method="BH"),
    cor_type = ifelse(R[ut] > 0, "Positive", "Negative"),
    stringsAsFactors = FALSE
  )
  
  sig <- if (allow_negative) subset(cor_df, weight >= rho_thr & q < q_thr)
  else                 subset(cor_df, rho    >= rho_thr & q < q_thr)
  
  if (nrow(sig) == 0) stop("筛选后无边；请放宽 rho_thr 或 q_thr。")

  g <- igraph::graph_from_data_frame(sig, directed=FALSE)
  g <- delete_vertices(g, V(g)[degree(g) == 0])

  if (!is.null(tax_input)) {
    tax0 <- if (is.character(tax_input)) {
      if (grepl("\\.csv$", tax_input, ignore.case=TRUE)) {
        read.csv(tax_input, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
      } else {
        read.table(tax_input, header=TRUE, sep="\t", check.names=FALSE,
                   stringsAsFactors=FALSE, quote="", comment.char="", fill=TRUE)
      }
    } else as.data.frame(tax_input, check.names=FALSE)
    
    V(g)$name <- .trim(V(g)$name)
    
    if (!is.null(rownames(tax0)) && !all(rownames(tax0) == as.character(seq_len(nrow(tax0))))) {
      tax0$.__rowid__ <- .trim(rownames(tax0))
    }
    
    id_candidates <- c("OTUID","ASV","FeatureID","ID","#OTU ID")
    idcol <- id_candidates[id_candidates %in% colnames(tax0)]
    if (length(idcol) == 0) idcol <- if (!is.null(tax0$.__rowid__)) ".__rowid__" else colnames(tax0)[1]
    
    tax0[[idcol[1]]] <- .trim(tax0[[idcol[1]]])
    tax0 <- tax0[!duplicated(tax0[[idcol[1]]]), , drop=FALSE]
    rownames(tax0) <- tax0[[idcol[1]]]
    
    idx <- match(V(g)$name, rownames(tax0))
    message("taxonomy 匹配到的节点数: ", sum(!is.na(idx)), " / ", vcount(g))
    
    tax_cols <- intersect(colnames(tax0), c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
    if (length(tax_cols) > 0 && any(!is.na(idx))) {
      for (nm in tax_cols) {
        vals <- rep(NA_character_, vcount(g))
        ok <- !is.na(idx)
        vals[ok] <- as.character(tax0[[nm]][idx[ok]])
        g <- igraph::set_vertex_attr(g, nm, value = vals)
      }
    }
  }

  cl <- cluster_louvain(g)
  V(g)$module <- as.character(membership(cl))
  deg_total <- degree(g)
  
  modules <- V(g)$module; nodes <- V(g)$name
  mods_unique <- unique(modules)
  adj_list <- adjacent_vertices(g, V(g))
  
  k_is_self <- numeric(vcount(g))
  k_is_by_mod <- vector("list", length=vcount(g))
  for (i in seq_along(nodes)) {
    nbrs <- names(adj_list[[i]])
    nbr_mods <- modules[match(nbrs, nodes)]
    tab <- table(nbr_mods)
    ki_vec <- setNames(numeric(length(mods_unique)), mods_unique)
    ki_vec[names(tab)] <- as.numeric(tab)
    k_is_by_mod[[i]] <- ki_vec
    k_is_self[i] <- ki_vec[modules[i]]
  }
  
  kstats <- lapply(mods_unique, function(m){
    idxm <- which(modules == m); kis <- k_is_self[idxm]
    c(mean = mean(kis), sd = stats::sd(kis))
  }); names(kstats) <- mods_unique
  
  Zi <- sapply(seq_along(nodes), function(i){
    m <- modules[i]; mu <- kstats[[m]]["mean"]; sdv <- kstats[[m]]["sd"]
    ifelse(is.na(sdv) || sdv == 0, 0, (k_is_self[i] - mu)/sdv)
  })
  Pi <- sapply(seq_along(nodes), function(i){
    ki <- deg_total[i]; if (ki == 0) 0 else 1 - sum((k_is_by_mod[[i]]/ki)^2, na.rm=TRUE)
  })
  
  V(g)$Zi <- Zi; V(g)$Pi <- Pi
  V(g)$role <- ifelse(Zi >= 2.5 & Pi < 0.62, "module_hub",
                      ifelse(Zi < 2.5 & Pi >= 0.62, "connector",
                             ifelse(Zi >= 2.5 & Pi >= 0.62, "network_hub","peripheral")))

  pos_num <- sum(E(g)$rho > 0); neg_num <- sum(E(g)$rho < 0)
  stats <- data.frame(
    nodes_num = vcount(g), edges_num = ecount(g),
    positive.corr_num  = pos_num, positive.corr_prop = pos_num/ecount(g),
    negative.corr_num  = neg_num, negative.corr_prop = neg_num/ecount(g),
    average_degree = mean(degree(g)),
    average_shortest_path_length = suppressWarnings(mean_distance(g, directed=FALSE)),
    network_diameter = suppressWarnings(diameter(g, directed=FALSE)),
    network_density  = edge_density(g),
    clustering_coefficient = transitivity(g, type="global"),
    modularity = modularity(g, membership(cl), directed=FALSE)
  )
  
  dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  f_nodes <- file.path(output_dir, sprintf("%s_nodes_%s.csv",  prefix, ts))
  f_edges <- file.path(output_dir, sprintf("%s_edges_%s.csv",  prefix, ts))
  f_stats <- file.path(output_dir, sprintf("%s_stats_%s.csv",  prefix, ts))
  f_graph <- file.path(output_dir, sprintf("%s_graph_%s.graphml", prefix, ts))
  
  nodes_df <- data.frame(
    id=V(g)$name, module=V(g)$module, degree=deg_total,
    Zi=V(g)$Zi, Pi=V(g)$Pi, role=V(g)$role, stringsAsFactors=FALSE
  )
  keep_attr <- setdiff(vertex_attr_names(g), c("name","module","Zi","Pi","role"))
  for (att in keep_attr) nodes_df[[att]] <- igraph::get.vertex.attribute(g, att)
  
  edges_df <- igraph::as_data_frame(g, what="edges")
  colnames(edges_df)[1:2] <- c("source","target")
  
  write.csv(nodes_df, f_nodes, row.names=FALSE)
  write.csv(edges_df, f_edges, row.names=FALSE)
  write.csv(stats,    f_stats, row.names=FALSE)
  igraph::write_graph(g, f_graph, format="graphml")
  
  list(graph=g, nodes=nodes_df, edges=edges_df, stats=stats,
       files=list(nodes=f_nodes, edges=f_edges, stats=f_stats, graphml=f_graph))
}

# Sample-specific subnetwork extraction + metrics
calc_graph_metrics <- function(gx){
  if (vcount(gx) < 2 || ecount(gx) < 1) {
    return(data.frame(
      nodes_num = vcount(gx),
      edges_num = ecount(gx),
      density = NA_real_,
      avg_degree = NA_real_,
      clustering = NA_real_,
      avg_path = NA_real_,
      diameter = NA_real_
    ))
  }
  data.frame(
    nodes_num = vcount(gx),
    edges_num = ecount(gx),
    density = edge_density(gx),
    avg_degree = mean(degree(gx)),
    clustering = transitivity(gx, type="global"),
    avg_path = suppressWarnings(mean_distance(gx, directed=FALSE)),
    diameter = suppressWarnings(diameter(gx, directed=FALSE))
  )
}

subgraph_by_sample <- function(g_global, otu_mat_sample_by_asv, sample_id, presence_thr = 0){
  if (!sample_id %in% rownames(otu_mat_sample_by_asv)) stop("sample not found: ", sample_id)
  present_nodes <- names(which(otu_mat_sample_by_asv[sample_id, ] > presence_thr))
  present_nodes <- intersect(present_nodes, V(g_global)$name)
  induced_subgraph(g_global, vids = V(g_global)[name %in% present_nodes])
}

# PARALLEL version
get_sample_network_metrics <- function(g_global, otu_mat_sample_by_asv,
                                       presence_thr = 0,
                                       export_graphml = FALSE,
                                       graphml_dir = "sample_subgraphs",
                                       workers = max(1, parallel::detectCores() - 1)) {
  if (export_graphml) dir.create(graphml_dir, showWarnings = FALSE, recursive = TRUE)
  
  old_plan <- future::plan()
  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(old_plan), add = TRUE)
  
  sids <- rownames(otu_mat_sample_by_asv)
  
  out <- future.apply::future_lapply(sids, function(sid){
    gs <- subgraph_by_sample(g_global, otu_mat_sample_by_asv, sid, presence_thr = presence_thr)
    m  <- calc_graph_metrics(gs)
    m$SampleID <- sid
    
    if (export_graphml) {
      f <- file.path(graphml_dir, paste0("subgraph_", sid, ".graphml"))
      igraph::write_graph(gs, f, format = "graphml")
    }
    m
  }, future.seed = TRUE)
  
  dplyr::bind_rows(out) %>% dplyr::select(SampleID, dplyr::everything())
}

# RUN

otu_file <- "Bacteria_otutab_rare.txt" 
tax_file <- "Bacteria_taxonomy.txt"       

first <- readLines(otu_file, n = 1)
is_header <- grepl("^#?OTU", first, ignore.case = TRUE) || grepl("Bacteria_ASV", first)

otu_raw <- read.table(
  otu_file,
  header = is_header,  
  row.names = 1,
  sep = "\t",         
  quote = "",
  comment.char = "",  
  check.names = FALSE
)
otu_mat <- t(otu_raw)  

res <- build_network(
  otu_input = otu_file,
  tax_input = tax_file,
  pre_filter = TRUE,
  prev_min = 0.20,
  mean_rel_min = 1e-4,
  pool_global = FALSE,
  rho_thr = 0.6,
  q_thr = 0.01,
  allow_negative = TRUE,
  output_dir = "Bacteria_Subnetwork_out",
  prefix = "Bacteria_Subnetwork"
)

g_global <- res$graph

overlap <- intersect(colnames(otu_mat), V(g_global)$name)
cat("Overlap between OTU columns and network nodes:", length(overlap), "\n")
if (length(overlap) < 10) stop("重叠太少：检查 V(g)$name 与 OTU 列名是否一致。")

otu_aligned <- otu_mat[, overlap, drop=FALSE]

metrics_df <- get_sample_network_metrics(
  g_global,
  otu_aligned,
  presence_thr = 0,
  export_graphml = FALSE,
  graphml_dir = "sample_subgraphs",
  workers = 24   
)

write.csv(metrics_df, "Bacteria_sample_subnetwork_metrics.csv", row.names=FALSE)

print(res$stats)
print(head(metrics_df))






library(ggplot2)

## Read and merge data
metadata <- read.delim("Bacteria_metadata_raw.txt", header = TRUE)

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
write.table(stat_res, "Bacteria_Subnetwork_Wilcoxon_BH.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## Color
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
  out_name <- paste0("Bacteria_Subnetwork_", v, "_boxplot.pdf")
  
  # Call the plotting function
  plot_metrics_df(
    data = net_meta,
    yvar = v,
    stat_table = stat_res,
    ylab = v,           # Directly use the network metric name as the Y-axis label
    out_pdf = out_name
  )
})


################################################################################################### Mean and SD Calculation
library(dplyr)

df <- read.csv("Bacteria_sample_subnetwork_metrics.csv", stringsAsFactors = FALSE)
metadata <- read.delim("Bacteria_metadata_raw.txt", stringsAsFactors = FALSE)

df2 <- merge(df, metadata, by = "SampleID")

metric_cols <- c("nodes_num", "edges_num", "density",
                 "avg_degree", "clustering", "avg_path", "diameter")

calc_mean <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  if (length(x) == 1) return(x)
  mean(x)
}

calc_sd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) <= 1) return(NA)
  sd(x)
}

result <- df2 %>%
  group_by(Group_Study) %>%
  summarise(across(all_of(metric_cols),
                   list(mean = calc_mean,
                        sd = calc_sd),
                   .names = "{.col}_{.fn}"),
            .groups = "drop")

write.table(result, "network_metrics_group_mean_sd.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)







































































