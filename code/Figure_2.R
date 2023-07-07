### Code for figure 2
### Raphael Eisenhofer 4/2023

##########################################################################################
## Load libraries, import data

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


#Subset by mouse and get MAG stats
MAG_stats %>%
  filter(str_detect(genome, "C04M3")) %>%
  group_by(seq_strat) %>%
  summarise(mean_n_contig = mean(n_contigs), mean_contam = mean(contamination),
            mean_comp = mean(completeness), n_MAGs = n(),
            n_90comp_MAGs = n_distinct(genome[completeness >= 90]))

MAG_stats %>%
  filter(str_detect(genome, "C13F5")) %>%
  group_by(seq_strat) %>%
  summarise(mean_n_contig = mean(n_contigs), mean_contam = mean(contamination),
            mean_comp = mean(completeness), n_MAGs = n(), 
            n_90comp_MAGs = n_distinct(genome[completeness >= 90]))


#Set theme and colours
theme_RE <- theme(
  text = element_text(family = "Arial"),
  axis.title.x = element_markdown(size = 15, family = "Arial"),
  axis.title.y = element_markdown(size = 15, family = "Arial"),
  axis.text.x = element_markdown(size = 12, family = "Arial"),
  axis.text.y = element_markdown(size = 12, family = "Arial"),
  panel.grid.major.x = element_blank(),
  panel.grid.minor = element_blank(),
  legend.text = element_text(size = 12),
  legend.title = element_blank())

colours <- c("ShortRead_20Gbp" = "#74a9cf", "ShortRead_40Gbp" = "#0570b0",
             "Hybrid_40Gbp" = "#fc8d59", "LongRead_20Gbp" ="#969696")


#Figure 2a: scatter plot
fig2a <- ggplot(MAG_stats, aes(completeness, contamination, colour = seq_strat, label = seq_strat)) +
  geom_point(aes(size=N50k), alpha = .6) +
  geom_smooth(se = FALSE) +
  scale_y_reverse() +
  scale_size(breaks = c(4000, 1000, 500, 100)) +
  theme_classic() +
  theme_RE +
  theme(legend.position = c(0.1, 0.15),
        legend.title = element_text(),
        legend.margin = margin(0, 0, 0, 0)
  ) +
  scale_color_manual(values = colours, name = "") +
  guides(colour = guide_legend(override.aes = list(linetype = c(0),
                                                   shape = c(16),
                                                   size = c(5))
  ),
  size = guide_legend(override.aes = list(colour = "grey")
  )
  ) +
  labs(x = "CheckM completeness", y = "CheckM contamination")



## N.B., to access the new column created in the pipe, enclose ggplot call
## with {} and reference it as '.' in the ggplot call.
fig2b <- MAG_stats %>%
  group_by(seq_strat) %>%
  arrange(n_contigs) %>%
  mutate(genome_number = 1:n()) %>%
  ungroup() %>%
  {
    ggplot(., aes(x = genome_number, y = n_contigs, fill = seq_strat)) +
      #  geom_col(aes(alpha = 0.6)) +
      geom_point(aes(alpha = 0.6, colour = seq_strat)) +
      geom_line(aes(colour = seq_strat)) +
      #  scale_y_continuous(breaks = seq(0, 900, 100)) +
      scale_x_continuous(limits = c(0, max(.$genome_number)+1), expand = c(0,0)) +
      scale_fill_manual(values = colours) +
      scale_colour_manual(values = colours) + 
      theme_minimal() +
      theme_RE +
      theme(
        #    panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_markdown(size = 10, family = "Arial"),
        legend.position = 0
      ) +
      labs(x = "Number of MAGs", y = "Number of contigs per MAG")
  }


fig2c <- MAG_stats %>%
  group_by(seq_strat) %>%
  arrange(n_contigs) %>%
  mutate(genome_number = 1:n()) %>%
  ungroup() %>%
  filter(genome_number <= 50) %>%
  {
    ggplot(., aes(x = genome_number, y = n_contigs, fill = seq_strat)) +
      #  geom_col(aes(alpha = 0.6)) +
      geom_point(aes(alpha = 0.6, colour = seq_strat)) +
      geom_line(aes(colour = seq_strat)) +
      #  scale_y_continuous(breaks = seq(0, 900, 100)) +
      scale_x_continuous(limits = c(0, 51), expand = c(0,0)) +
      scale_y_continuous(breaks = c(0, 10, seq(5, 100, 10))) +
      scale_fill_manual(values = colours) +
      scale_colour_manual(values = colours) + 
      theme_minimal() +
      theme_RE +
      theme(
        #    panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_markdown(size = 10, family = "Arial"),
        legend.position = 0
      ) +
      labs(x = "Number of MAGs", y = "Number of contigs per MAG")
  }


#Stitch together with patchwork
fig2a + (fig2b / fig2c) +
  plot_layout(widths = c(10, 5)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_markdown(size = 18, face = "bold")) 

ggsave(filename = "figures/Figure_2.png", width = 15, height = 10, units = "in")
