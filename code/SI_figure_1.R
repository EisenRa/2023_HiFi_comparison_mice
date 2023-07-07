################################################################################
### SI Figure 1: Raincloud plots of CheckM completeness and contamination
### Raphael Eisenhofer 4/2023

library(tidyverse)
library(systemfonts)
library(ggtext)
library(ggforce)
library(ggdist)
library(ggpubr)
library(patchwork)

#Get MAGs found in all sequencing strategies
n_contigs <- read_delim("3_Outputs/mag_stats/n_contigs_per_bin.tsv", col_names = c("genome", "n_contigs")) %>%
  mutate(genome = str_replace_all(genome, ".gz", ""))
genome_info <- read_delim("3_Outputs/mag_stats/genomeInformation.csv") %>%
  mutate(genome = str_replace_all(genome, ".gz", ""))

#Merge data into a single frame, mutate new columns
MAG_stats <- left_join(n_contigs, genome_info, by = "genome") %>%
  mutate_at(., "genome", str_replace, ".fa", "") %>%
  mutate(N50k = N50 / 1000,
         seq_strat = case_when(
           str_detect(genome, "hifiasm") ~ "LongRead_20Gbp",
           str_detect(genome, "subsampled_20Gb_Hybrid") ~ "Hybrid_40Gbp",
           str_detect(genome, "subsampled_20Gb") ~ "ShortRead_20Gbp",
           str_detect(genome, "subsampled_40Gb") ~ "ShortRead_40Gbp")) %>%
  mutate(seq_strat = fct_relevel(seq_strat, c("ShortRead_20Gbp", "ShortRead_40Gbp", 
                                              "Hybrid_40Gbp", "LongRead_20Gbp")))


break_labels <- c("ShortRead<br>20Gbp", "ShortRead<br>40Gbp", 
                  "Hybrid<br>40Gbp", "LongRead<br>20Gbp")

si_fig1a <- ggplot(MAG_stats, aes(seq_strat, completeness, 
                                  colour = seq_strat, label = seq_strat)) +
  stat_halfeye(
    adjust = .5,
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA
  ) + 
  geom_point(
    aes(colour = seq_strat),
    size = 4,
    alpha = .6,
    position = position_jitter(seed = 1, width = .1)
  ) + 
  stat_summary(
    fun = "mean", 
    geom = "crossbar", 
    colour = "red", 
    width = 0.2
  ) +
  stat_summary(
    fun = "mean",
    geom = "text",
    aes(label = round(..y.., 2)),
    hjust = 1.75,
    colour = "red"
  ) +  
  theme_light() +
  theme_RE +
  theme(
    #        legend.position = c(0.85, 0.15),
    legend.box.background = element_rect(size = 0.5),
    legend.margin = margin(-5, 5, 0, 0),
    axis.title.x = element_blank(),
    legend.position = "none") +
  scale_color_manual(values = colours) +
  scale_x_discrete(labels = break_labels) +
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  labs(x = "<br>Sequencing strategy", 
       y = "CheckM completeness"
  )


si_fig1b <- ggplot(MAG_stats, aes(seq_strat, contamination, 
                                  colour = seq_strat, label = seq_strat)) +
  stat_halfeye(
    adjust = .5,
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA
  ) + 
  geom_point(
    aes(colour = seq_strat),
    size = 4,
    alpha = .6,
    position = position_jitter(seed = 1, width = .1)
  ) + 
  stat_summary(
    fun = "mean", 
    geom = "crossbar", 
    colour = "red", 
    width = 0.2
  ) +
  stat_summary(
    fun = "mean",
    geom = "text",
    aes(label = round(..y.., 2)),
    hjust = 2,
    colour = "red"
  ) +
  theme_light() +
  theme_RE +
  theme(
    #        legend.position = c(0.35, 0.875),
    legend.box.background = element_rect(size = 0.5),
    legend.margin = margin(-5, 5, 0, 0),
    axis.title.x = element_blank(),
    legend.position = "none") +
  scale_color_manual(values = colours) +
  scale_x_discrete(labels = break_labels) +
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  labs(x = "<br>Sequencing strategy", 
       y = "CheckM contamination"
  )


#Stitch together with patchwork
(si_fig1a | si_fig1b) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_markdown(size = 18, face = "bold")) 

ggsave(filename = "figures/SI_figure_1.png", width = 10, height = 5, units = "in")
