

library(dplyr)
library(ggplot2)
library(maps)
library(tidyverse)
library(readxl)

world.dat <- map_data("world")

ggplot() +
  geom_polygon(data = world.dat, aes(x = long, y = lat, group = group),
               fill = "#d4d4d4", color = NA) + 
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(add = c(0, 0))) +
  theme(
    panel.background = element_rect(fill = "#e3e9ef", color = NA),  
    plot.background = element_rect(fill = "#e3e9ef", color = NA),   
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank()
  ) +
  labs(x = NULL, y = NULL) -> world.map


df <- read.csv(file.choose())

p <- world.map +
  geom_point(data = df,
             aes(x = Longitude,
                 y = Latitude,
                 size = Size,
                 fill = Type),
             color = "black",
             shape = 21,
             alpha = 0.5,
             stroke = 0.3) +
  scale_fill_manual(
    values = c("Bacteria" = "#448DCD", "Fungi" = "#E37B6D"),
    name = "Type"
  ) +
  scale_size_continuous(name = "number", range = c(2, 10)) +
  theme(
    legend.position = c(0.12, 0.42)
  )
p


ggsave("map.pdf",
       p,
       width = 8,
       height = 4)

# 8*4

