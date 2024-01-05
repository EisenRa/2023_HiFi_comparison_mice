library(tidyverse)
library(ggpattern)

rrnas <- read_delim("data/mag_stats/barrnap_results.tsv", col_names = c("length", "mag"), delim = '\t') %>%
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
  ) %>%
  pivot_wider(id_cols = mag, names_from = type, values_from = type, values_fn = list(type = length))

trnas <- read_delim("data/mag_stats/trna_stats.tsv", col_names = c("mag", "num_trnas"), delim = '\t')

genome_info <- read_delim("data/mag_stats/genomeInformation.csv") %>%
  mutate(genome = str_replace_all(genome, ".fa.gz", "")) %>%
  rename(mag = genome)

df <- genome_info %>%
  left_join(., trnas, by = join_by("mag")) %>%
  left_join(., rrnas, by = join_by("mag")) %>%
  mutate(seq_strat = case_when(
    str_detect(mag, "hifiasm") ~ "LongRead_20Gbp",
    str_detect(mag, "subsampled_20Gb_Hybrid") ~ "Hybrid_40Gbp",
    str_detect(mag, "subsampled_20Gb") ~ "ShortRead_20Gbp",
    str_detect(mag, "subsampled_40Gb") ~ "ShortRead_40Gbp")) %>%
  mutate(seq_strat = fct_relevel(seq_strat, c("ShortRead_20Gbp", "ShortRead_40Gbp", 
                                                "Hybrid_40Gbp", "LongRead_20Gbp"))) %>%
  mutate(mimag = case_when(
    completeness >= 90 & contamination < 5 & `16S` >= 1 & `5S` >= 1 & `23S` >= 1 & num_trnas >=18 ~ "high-quality",
    .default = "medium-quality"
  ))


# Plotting
colours <- c("short_20gbp" = "#74a9cf", "short_40gbp" = "#0570b0",
             "hybrid" = "#fc8d59", "long_read" ="#969696")

df %>%
  filter(mimag == "high-quality") %>%
  summarise(n = n(), .by = seq_strat)

