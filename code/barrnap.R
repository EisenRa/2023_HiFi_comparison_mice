library(tidyverse)
library(ggpattern)

df <- read_delim("data/mag_stats/barrnap_results.tsv", col_names = c("length", "mag"), delim = '\t') %>%
  mutate(length = str_replace_all(length, " \\(partial\\)", ""),
         type = word(length, -3),
         length = as.numeric(str_extract(length, "(?<=L=)\\d+")),
         assembly = case_when(str_detect(mag, "hifiasm") ~ "long_read",
                              str_detect(mag, "Hybrid") ~ "hybrid",
                              str_detect(mag, "20Gb_bin") ~ "short_20gbp",
                              str_detect(mag, "40Gb_bin") ~ "short_40gbp"
                              ),
         length_class = case_when(type == "16S" & length > 1450 ~ "full",
                                  type == "16S" & length < 1450 ~ "incomplete")
         )


colours <- c("short_20gbp" = "#74a9cf", "short_40gbp" = "#0570b0",
             "hybrid" = "#fc8d59", "long_read" ="#969696")

df %>%
  filter(type == "16S") %>%
  ggplot(aes(x = factor(assembly, level = c("long_read", "hybrid", "short_20gbp", "short_40gbp")), fill = assembly)) +
  geom_bar_pattern(mapping = aes(pattern = length_class), pattern_density = 1, pattern_fill = "black", pattern_fill2 = "black") +
  scale_fill_manual(values = colours, name = "") +
  scale_pattern_manual(values = c("full" = "none", "incomplete" = "stripe")) +
  theme_classic() + 
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(x = "Assembly type", y = "Number of 16S rRNA genes recovered")

ggsave("figures/SI_figure_3.png", width = 5, height = 5)  
