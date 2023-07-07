library(tidyverse)
library(janitor)

c04m3_long_all_checkm <- read_delim("HIFI_ASM/C04M3_longread_hifiasm_ALL.fa_checkm.tsv")
c04m3_long_circ_checkm <- read_delim("HIFI_ASM/C04M3_longread_hifiasm_CIRCULAR.fa_checkm.tsv")
c13m5_long_all_checkm <- read_delim("HIFI_ASM/C13M5_longread_hifiasm_ALL.fa_checkm.tsv")
c13m5_long_circ_checkm <- read_delim("HIFI_ASM/C13M5_longread_hifiasm_CIRCULAR.fa_checkm.tsv")

contig_stats <- read_delim("HIFI_ASM/hifiasm_meta_contig_stats.tsv", col_names = c("bin_id", "ncontigs", "length")) %>%
  mutate(bin_id = str_replace_all(bin_id, ".fa", ""))

cdb <- read_delim("HIFI_ASM/Cdb.csv") %>%
  rename("bin_id" = genome) %>%
  mutate(bin_id = str_replace_all(bin_id, ".fa", ""))
  select(bin_id, secondary_cluster)

widb <- read_delim("HIFI_ASM/dRep98_Widb.csv") %>%
  rename("bin_id" = genome) %>%
  mutate(bin_id = str_replace_all(bin_id, ".fa", "")) %>%
  select(bin_id, score, cluster, cluster_members, furthest_cluster_member)

df <- rbind(c04m3_long_all_checkm, c04m3_long_circ_checkm, c13m5_long_all_checkm, c13m5_long_circ_checkm) %>%
  clean_names() %>%
  select(bin_id, number_markers, completeness, contamination, strain_heterogeneity) %>%
  inner_join(., contig_stats, by = "bin_id")

bact_contigs <- cdb %>%
  inner_join(., df, by = "bin_id")

winners <- widb %>%
  inner_join(., df, by = "bin_id")
