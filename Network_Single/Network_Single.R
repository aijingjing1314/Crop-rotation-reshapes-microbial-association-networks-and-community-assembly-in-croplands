

library(ggplot2)
library(dplyr)

Network_single <- read.table("Network_Single_Only8.txt", header = TRUE, fill = TRUE, stringsAsFactors = FALSE)


# Claculate Mean and SE
summary_data <- Network_single %>%
  summarise(
    Bacteria_mean_nodes = mean(Bacteria_nodes_num_change, na.rm = TRUE),
    Bacteria_se_nodes = sd(Bacteria_nodes_num_change, na.rm = TRUE) / sqrt(sum(!is.na(Bacteria_nodes_num_change))),
    Bacteria_mean_edges = mean(Bacteria_edges_num_change, na.rm = TRUE),
    Bacteria_se_edges = sd(Bacteria_edges_num_change, na.rm = TRUE) / sqrt(sum(!is.na(Bacteria_edges_num_change))),
    Fungi_mean_nodes = mean(Fungi_nodes_num_change, na.rm = TRUE),
    Fungi_se_nodes = sd(Fungi_nodes_num_change, na.rm = TRUE) / sqrt(sum(!is.na(Fungi_nodes_num_change))),
    Fungi_mean_edges = mean(Fungi_edges_num_change, na.rm = TRUE),
    Fungi_se_edges = sd(Fungi_edges_num_change, na.rm = TRUE) / sqrt(sum(!is.na(Fungi_edges_num_change)))
  )

plot_Bacteria_nodes <- data.frame(
  Metric = "Bacteria_Nodes Number Change",
  Mean = summary_data$Bacteria_mean_nodes,
  SE = summary_data$Bacteria_se_nodes
)

plot_Bacteria_edges <- data.frame(
  Metric = "Bacteria_Edges Number Change",
  Mean = summary_data$Bacteria_mean_edges,
  SE = summary_data$Bacteria_se_edges
)

plot_Fungi_nodes <- data.frame(
  Metric = "Fungi_Nodes Number Change",
  Mean = summary_data$Fungi_mean_nodes,
  SE = summary_data$Fungi_se_nodes
)

plot_Fungi_edges <- data.frame(
  Metric = "Fungi_Edges Number Change",
  Mean = summary_data$Fungi_mean_edges,
  SE = summary_data$Fungi_se_edges
)


p1 <- ggplot(plot_Bacteria_nodes, aes(x = Metric, y = Mean)) +
  geom_col(fill = "grey", width = 0.6, color = "black") + 
  geom_errorbar(aes(ymin = Mean, ymax = Mean + SE),
                width = 0.15, linewidth = 0.8) +
  theme_classic() +
  labs(title = "Bacteria_Nodes Number Change (Mean ± SE)", x = NULL, y = "Change (%)") +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    aspect.ratio = 1.5 
  )
print(p1)
ggsave(filename = "Bacteria_Nodes_Change.pdf", plot = p1, width = 3, height = 5)

p2 <- ggplot(plot_Bacteria_edges, aes(x = Metric, y = Mean)) +
  geom_col(fill = "grey", width = 0.6, color = "black") + 
  geom_errorbar(aes(ymin = Mean, ymax = Mean + SE),
                width = 0.15, linewidth = 0.8) +
  theme_classic() +
  labs(title = "Bacteria_Nodes Number Change (Mean ± SE)", x = NULL, y = "Change (%)") +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    aspect.ratio = 1.5 
  )
print(p2)
ggsave(filename = "Bacteria_Edges_Change.pdf", plot = p2, width = 3, height = 5)

p3 <- ggplot(plot_Fungi_nodes, aes(x = Metric, y = Mean)) +
  geom_col(fill = "grey", width = 0.6, color = "black") + 
  geom_errorbar(aes(ymin = Mean, ymax = Mean + SE),
                width = 0.15, linewidth = 0.8) +
  theme_classic() +
  labs(title = "Fungi_Nodes Number Change (Mean ± SE)", x = NULL, y = "Change (%)") +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    aspect.ratio = 1.5 
  )
print(p3)
ggsave(filename = "Fungi_Nodes_Change.pdf", plot = p3, width = 3, height = 5)

p4 <- ggplot(plot_Fungi_edges, aes(x = Metric, y = Mean)) +
  geom_col(fill = "grey", width = 0.6, color = "black") + 
  geom_errorbar(aes(ymin = Mean, ymax = Mean + SE),
                width = 0.15, linewidth = 0.8) +
  theme_classic() +
  labs(title = "Fungi_Nodes Number Change (Mean ± SE)", x = NULL, y = "Change (%)") +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    aspect.ratio = 1.5 
  )
print(p4)
ggsave(filename = "Fungi_Edges_Change.pdf", plot = p4, width = 3, height = 5)






####################################################################  sensitivity analysis

library(ggplot2)
library(dplyr)

Network_single <- read.table("Network_Single_Only8_sensitivity.txt", header = TRUE, fill = TRUE, stringsAsFactors = FALSE)


# Claculate Mean and SE
summary_data <- Network_single %>%
  summarise(
    Bacteria_mean_nodes = mean(Bacteria_nodes_num_change, na.rm = TRUE),
    Bacteria_se_nodes = sd(Bacteria_nodes_num_change, na.rm = TRUE) / sqrt(sum(!is.na(Bacteria_nodes_num_change))),
    Bacteria_mean_edges = mean(Bacteria_edges_num_change, na.rm = TRUE),
    Bacteria_se_edges = sd(Bacteria_edges_num_change, na.rm = TRUE) / sqrt(sum(!is.na(Bacteria_edges_num_change))),
    Fungi_mean_nodes = mean(Fungi_nodes_num_change, na.rm = TRUE),
    Fungi_se_nodes = sd(Fungi_nodes_num_change, na.rm = TRUE) / sqrt(sum(!is.na(Fungi_nodes_num_change))),
    Fungi_mean_edges = mean(Fungi_edges_num_change, na.rm = TRUE),
    Fungi_se_edges = sd(Fungi_edges_num_change, na.rm = TRUE) / sqrt(sum(!is.na(Fungi_edges_num_change)))
  )

plot_Bacteria_nodes <- data.frame(
  Metric = "Bacteria_Nodes Number Change",
  Mean = summary_data$Bacteria_mean_nodes,
  SE = summary_data$Bacteria_se_nodes
)

plot_Bacteria_edges <- data.frame(
  Metric = "Bacteria_Edges Number Change",
  Mean = summary_data$Bacteria_mean_edges,
  SE = summary_data$Bacteria_se_edges
)

plot_Fungi_nodes <- data.frame(
  Metric = "Fungi_Nodes Number Change",
  Mean = summary_data$Fungi_mean_nodes,
  SE = summary_data$Fungi_se_nodes
)

plot_Fungi_edges <- data.frame(
  Metric = "Fungi_Edges Number Change",
  Mean = summary_data$Fungi_mean_edges,
  SE = summary_data$Fungi_se_edges
)


p1 <- ggplot(plot_Bacteria_nodes, aes(x = Metric, y = Mean)) +
  geom_col(fill = "grey", width = 0.6, color = "black") + 
  geom_errorbar(aes(ymin = Mean, ymax = Mean + SE),
                width = 0.15, linewidth = 0.8) +
  theme_classic() +
  labs(title = "Bacteria_Nodes Number Change (Mean ± SE)", x = NULL, y = "Change (%)") +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    aspect.ratio = 1.5 
  )
print(p1)
ggsave(filename = "Bacteria_Nodes_Change_sensitivity.pdf", plot = p1, width = 3, height = 5)

p2 <- ggplot(plot_Bacteria_edges, aes(x = Metric, y = Mean)) +
  geom_col(fill = "grey", width = 0.6, color = "black") + 
  geom_errorbar(aes(ymin = Mean, ymax = Mean + SE),
                width = 0.15, linewidth = 0.8) +
  theme_classic() +
  labs(title = "Bacteria_Nodes Number Change (Mean ± SE)", x = NULL, y = "Change (%)") +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    aspect.ratio = 1.5 
  )
print(p2)
ggsave(filename = "Bacteria_Edges_Change_sensitivity.pdf", plot = p2, width = 3, height = 5)

p3 <- ggplot(plot_Fungi_nodes, aes(x = Metric, y = Mean)) +
  geom_col(fill = "grey", width = 0.6, color = "black") + 
  geom_errorbar(aes(ymin = Mean, ymax = Mean + SE),
                width = 0.15, linewidth = 0.8) +
  theme_classic() +
  labs(title = "Fungi_Nodes Number Change (Mean ± SE)", x = NULL, y = "Change (%)") +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    aspect.ratio = 1.5 
  )
print(p3)
ggsave(filename = "Fungi_Nodes_Change_sensitivity.pdf", plot = p3, width = 3, height = 5)

p4 <- ggplot(plot_Fungi_edges, aes(x = Metric, y = Mean)) +
  geom_col(fill = "grey", width = 0.6, color = "black") + 
  geom_errorbar(aes(ymin = Mean, ymax = Mean + SE),
                width = 0.15, linewidth = 0.8) +
  theme_classic() +
  labs(title = "Fungi_Nodes Number Change (Mean ± SE)", x = NULL, y = "Change (%)") +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    aspect.ratio = 1.5 
  )
print(p4)
ggsave(filename = "Fungi_Edges_Change_sensitivity.pdf", plot = p4, width = 3, height = 5)

