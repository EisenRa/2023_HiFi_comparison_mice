library(tidyverse)
library(janitor)

short40gb_bins <- 191
short20gb_bins <- 142
hybrid_bins <- 173
long_bins <- 113

num_bins <- tibble(
  seq_strat = c("short_40gbp", "short_20gbp", "hybrid", "long_read"),
  n_bins = c(short40gb_bins, short20gb_bins, hybrid_bins, long_bins)
)

colours <- c("short_20gbp" = "#74a9cf", "short_40gbp" = "#0570b0",
             "hybrid" = "#fc8d59", "long_read" ="#969696")

df <- read_delim("3_Outputs/mag_stats/antismash_summary.tsv", 
                 delim="\t", 
                 col_names = c("region", "product", "region_length", 
                               "n_genes", "contig_edge")) %>%
  mutate(product = str_replace_all(product, ".*/product=", ""),
         product = str_replace_all(product, "\"", ""),
         contig_edge = str_replace_all(contig_edge, ".*/contig_edge=", ""),
         contig_edge = str_replace_all(contig_edge, "\"", ""),
         seq_strat = case_when(str_detect(region, "hifiasm") ~ "long_read",
                               str_detect(region, "subsampled_20Gb_bin") ~ "short_20gbp",
                               str_detect(region, "subsampled_40Gb") ~ "short_40gbp",
                               str_detect(region, "Hybrid") ~ "hybrid"),
         mags_per_group = case_when(seq_strat == "long_read" ~ long_bins,
                                    seq_strat == "hybrid" ~ hybrid_bins,
                                    seq_strat == "short_20gbp" ~ short20gb_bins,
                                    seq_strat == "short_40gbp" ~ short40gb_bins)
  ) %>%
  separate_wider_delim(region, "_antismash/", names = c("mag", "region"))



num_bgcs <- df %>%
  group_by(seq_strat) %>%
  summarise(n_bgcs = n()) %>%
  inner_join(., num_bins, by = "seq_strat") %>%
  mutate(normalised = n_bgcs / n_bins)

region_length <- df %>%
  group_by(seq_strat) %>%
  summarise(bgcs_length = mean(region_length))

product <- df %>%
  group_by(seq_strat, product) %>%
  summarise(n = n())

product %>%
  ggplot(aes(y = n, x = product, fill = seq_strat)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = colours, name = "") +
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle = 90, size = 12, vjust = 0.2, hjust = 0.95),
    axis.text.y = element_text(size = 12),
    legend.position = c(0.1, 0.9)
    ) +
  labs(x = "biosynthetic gene cluster type", y = "number (n)")

ggsave("figures/figure_s2.png")

df %>%
  group_by(seq_strat, contig_edge) %>%
  summarise(n = n()) %>%
  inner_join(., num_bgcs, by = "seq_strat") %>%
  mutate(percent =  (n / n_bgcs) * 100) %>%
  inner_join(., region_length, by = "seq_strat")



