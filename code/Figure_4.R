### Code for figure 4
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
clusters <- read_delim("3_Outputs/mag_stats/Cdb.csv") %>%
  mutate(genome = str_replace_all(genome, ".gz", ""))
n_contigs <- read_delim("3_Outputs/mag_stats/n_contigs_per_bin.tsv", col_names = c("genome", "n_contigs")) %>%
  mutate(genome = str_replace_all(genome, ".gz", ""))
genome_info <- read_delim("3_Outputs/mag_stats/genomeInformation.csv") %>%
  mutate(genome = str_replace_all(genome, ".gz", ""))


MAG_stats <- left_join(clusters, n_contigs, by = "genome")
MAG_stats <- left_join(MAG_stats, genome_info, by = "genome") %>%
  select(-comparison_algorithm, -cluster_method, -threshold, -centrality) %>%
  mutate_at(., "genome", str_replace, ".fa", "") %>%
  mutate(seq_strat = case_when(
    str_detect(genome, "hifiasm") ~ "LongRead_20Gbp",
    str_detect(genome, "subsampled_20Gb_Hybrid") ~ "Hybrid_40Gbp",
    str_detect(genome, "subsampled_20Gb") ~ "ShortRead_20Gbp",
    str_detect(genome, "subsampled_40Gb") ~ "ShortRead_40Gbp")) %>%
  mutate(seq_strat = fct_relevel(seq_strat, c("ShortRead_20Gbp", "ShortRead_40Gbp", 
                                              "Hybrid_40Gbp", "LongRead_20Gbp")))


##########################################################################################
## Import count data, find primary clusters that don't have long-read derived assemblies,
## normalise and filter count table, then create a phyloseq object

library(phyloseq)
library(scales)
library(ggh4x)

#Load in the count table:
count_table <- read_delim("3_Outputs/mag_stats/coverm_mag_stats.tsv", delim = '\t')

#Get list of dereplicated MAGs, join to MAG_stats table
dereplicated_mags <- count_table %>%
  select(Genome) %>%
  mutate(dereplicated_mag = "TRUE")

#Join tables then create a new column to indicate if genome is long read-derived
MAG_stats <- left_join(MAG_stats, dereplicated_mags, by = c("genome" = "Genome")) %>%
  mutate(
    LongRead_MAG = case_when(
      str_detect(seq_strat, "LongRead") ~ "TRUE", TRUE ~ "FALSE"),
    Hybrid_MAG = case_when(
      str_detect(seq_strat, "Hybrid") ~ "TRUE", TRUE ~ "FALSE"),
    Short_read_20 = case_when(
      str_detect(seq_strat, "ShortRead_20Gbp") ~ "TRUE", TRUE ~ "FALSE"),
    Short_read_40 = case_when(
      str_detect(seq_strat, "ShortRead_40Gbp") ~"TRUE", TRUE ~ "FALSE")
  )



#Now, group by primary cluster, then keep only genomes from primary clusters that have
#a long read-derived MAG in them, and create a new column indicating this.
MAGs_with_long_read_assemblies_in_primary_cluster <- MAG_stats %>%
  group_by(primary_cluster) %>%
  filter(any(LongRead_MAG == "TRUE")) %>%
  ungroup() %>%
  select(genome, -primary_cluster) %>%
  mutate(LongRead_MAG_in_final_cluster = "TRUE")

MAGs_with_hybrid_assemblies_in_primary_cluster <- MAG_stats %>%
  group_by(primary_cluster) %>%
  filter(any(Hybrid_MAG == "TRUE")) %>%
  ungroup() %>%
  select(genome, -primary_cluster) %>%
  mutate(Hybrid_MAG_in_final_cluster = "TRUE")

MAGs_with_short_read_20_assemblies_in_primary_cluster <- MAG_stats %>%
  group_by(primary_cluster) %>%
  filter(any(Short_read_20 == "TRUE")) %>%
  ungroup() %>%
  select(genome, -primary_cluster) %>%
  mutate(ShortRead20_MAG_in_final_cluster = "TRUE")

MAGs_with_short_read_40_assemblies_in_primary_cluster <- MAG_stats %>%
  group_by(primary_cluster) %>%
  filter(any(Short_read_40 == "TRUE")) %>%
  ungroup() %>%
  select(genome, -primary_cluster) %>%
  mutate(ShortRead40_MAG_in_final_cluster = "TRUE")

#Join these tables, then replace NAs with 'FALSE'. 
MAG_stats_final <- left_join(MAG_stats, MAGs_with_long_read_assemblies_in_primary_cluster,
                             by = "genome") %>%
  left_join(., MAGs_with_hybrid_assemblies_in_primary_cluster, by = "genome") %>%
  left_join(., MAGs_with_short_read_20_assemblies_in_primary_cluster, by = "genome") %>%
  left_join(., MAGs_with_short_read_40_assemblies_in_primary_cluster, by = "genome") %>%
  replace(is.na(.), "FALSE")

#Create lists of dereplicated MAGs that do/do not have long read-derived MAGs in their 
#primary cluster
MAGs_that_long_reads_missed <- MAG_stats_final %>%
  filter(dereplicated_mag == "TRUE" & LongRead_MAG_in_final_cluster == "FALSE") %>%
  select(genome)

MAGs_that_long_reads_missed_vector <- MAGs_that_long_reads_missed$genome

MAGs_with_long_read_representative <- MAG_stats_final %>%
  filter(dereplicated_mag == "TRUE" & LongRead_MAG_in_final_cluster == "TRUE") %>%
  select(genome)

MAGs_with_long_read_representative_vector <- MAGs_with_long_read_representative$genome

#Normalise the count table by MAG size and filter by evenness of coverage (detection):
#MAGs need reads covering at least 30% of the genome to be detected.
norm_count_table <- count_table %>% 
  pivot_longer(!Genome,
               names_to = c("sample", ".value"),
               names_sep = " ") %>%
  dplyr::rename("read_count" = "Read") %>%
  dplyr::rename("covered_fraction" = "Covered") %>%
  dplyr::rename("MAG_length" = "Length") %>%
  mutate(read_count = if_else(covered_fraction < 0.3, 0, read_count)) %>%
  mutate(read_count = read_count*200/MAG_length) %>%
  select(Genome, sample, read_count) %>%
  pivot_wider(names_from = sample, values_from = read_count)

#Create a phyloseq object
taxonomy <- read.delim("3_Outputs/mag_stats/gtdbtk.bac120.summary.tsv", sep = '\t') %>%
  mutate(user_genome = str_replace_all(user_genome, "_renamed.fa", ""))

## Coerce gtdb-tk output to phyloseq taxonomy table
# Remove the "p__" etc. designations from the classification string
taxonomyclean <- taxonomy %>% 
  mutate(classification = str_replace_all(classification, ".__", "")) %>%
  # Split the classification string into columns for each taxonnomic rank
  separate(col = classification, sep = ";", into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
# Remove first column
taxonomyclean1 <- taxonomyclean[,-1]
# Set MAG name as row names of the dataframe
rownames(taxonomyclean1) <- taxonomyclean[,1]

# Get the necessary columns, convert to matrix
taxmatrix <- as.matrix(taxonomyclean1[1:7])

taxtable <- tax_table(taxmatrix) 

## Coerce count table into a phylseq OTU table
# Convert to dataframe
MAGdf <- as.data.frame(norm_count_table)
# Remove first column (MAG names)
MAGdfrows <- MAGdf[,-1]
# Set MAG names as row names of the dataframe
rownames(MAGdfrows) <- MAGdf[,1] 

MAGotu <- otu_table(MAGdfrows, taxa_are_rows = T)

## Coerce metadata into metadata object

#Finally, create the phyloseq object
ps <- phyloseq(MAGotu, taxtable)

##########################################################################################
## Create heatmap of 20 Mbp subsampled data mapped to the dereplicated MAG catalogue

#Convert counts to relative abundance
ps.ra <- transform_sample_counts(ps, function(x) x / sum(x) )

#What's the relative abundance of MAGs without a long-read assembly representative?
percent(mean(sample_sums(prune_taxa(MAGs_that_long_reads_missed_vector, ps.ra))))


#Add sample temporal information so we can order samples by time
mouse_sample_times <- read_delim("3_Outputs/mouse_sample_times.txt") %>%
  select(!date) %>%
  separate(., sample_time, into = c("mouse_id", "order_number"), sep = "_") %>%
  mutate(sample = str_replace_all(sample, " ", ""))



df <- psmelt(ps.ra) %>%
  inner_join(., MAG_stats_final, by = c("OTU" = "genome")) %>%
  mutate(
    sample_id = case_when(
      str_detect(Sample, "C04M3") ~ "C04M3",
      str_detect(Sample, "C13F5") ~ "C13F5"),
    Sample = str_replace_all(Sample, "_ss_dRepd_mags", "")
  ) 

#Create a new column with Phylum__family__genus
df <- df %>% unite(., col = phy_fam_gen, Phylum, Family, Genus, sep = "  |  ", remove = FALSE)
df <- df %>% unite(., col = fam_gen, Family, Genus, sep = "  |  ", remove = FALSE)


#Join em up chief
df <- df %>%
  left_join(., mouse_sample_times, by = c("Sample" = "sample")) %>%
  unite(., col = OTU, Phylum, OTU, sep = "_", remove = FALSE) %>%
  mutate(LongRead_MAG_in_final_cluster = str_replace_all(LongRead_MAG_in_final_cluster, 
                                                         "TRUE", 
                                                         "MAGs in all strategies")) %>%
  mutate(LongRead_MAG_in_final_cluster = str_replace_all(LongRead_MAG_in_final_cluster, 
                                                         "FALSE", 
                                                         "MAGs not found by long-reads"))


#HEAT MAP HEAT MAP!
fig3a <- df %>%
  ggplot(aes(x = reorder(order_number, as.numeric(order_number)), 
             y = reorder(OTU, Abundance), 
             fill = Abundance)) +
  geom_tile() +
  facet_nested(rows = vars(LongRead_MAG_in_final_cluster),
               cols = vars(sample_id),
               scales = "free",
               space = "free",
               switch = "both",
#               margins = c(0,0)
  ) +
  #Find a way to manually set gradient!  
  scale_fill_gradientn(colours = c("white", "grey", "black"),
                       values = c(0, 0.5, 1),
  ) +
  theme_classic() +
  theme_RE + 
  theme(
    axis.text.x = element_markdown(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    strip.text.y = element_blank(),
    strip.text.x = element_markdown(size = 16),
    legend.position = "left",
    legend.title = element_markdown()
    #    legend.margin = margin(0, 0, -40, 0)
  ) +
  labs(x = "Sample day")


## Add the additional tiles to the heatmap

theme_tiles <- theme(
  axis.text.y = element_blank(),
  axis.text.x = element_markdown(size = 0.1),
  axis.title.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_markdown(size = 12),
  strip.text = element_blank(),
  axis.line.x = element_blank(),
  legend.position = "none",
  axis.title = element_blank(),
  axis.ticks.length = unit(0, "pt")
)


fig3_taxa_labels <- df %>%
  ggplot(aes(x = "", y = reorder(OTU, Abundance), 
             fill = Phylum)) +
  geom_tile(colour = "black") +
  scale_fill_manual(values = c("Bacteroidota" = "#4575b4", "Verrucomicrobiota" = "#91bfdb",
                               "Firmicutes_A" = "#d73027", "Firmicutes" = "#fc8d59",
                               "Firmicutes_B" = "#fee090", "Proteobacteria" = "#ffffbf",
                               "Actinobacteriota" = "#e0f3f8")) +
  facet_nested(rows = vars(LongRead_MAG_in_final_cluster),
               scales = "free_y",
               space = "free_y",
               switch = "both") +
  scale_x_discrete(expand = c(0, 0), position = "left") +
  theme_classic() +
  theme_RE +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    strip.text.y = element_markdown(size = 14),
    legend.text = element_markdown(size = 10),
    legend.position = "left"
  )
theme_tiles +
  xlab("Phylum classification of MAGs")



fig3_shortread40 <- df %>%
  ggplot(aes(x = "", y = reorder(OTU, Abundance), 
             fill = ShortRead40_MAG_in_final_cluster)) +
  geom_tile(colour = "white") +
  scale_fill_manual(values=c("FALSE"="white","TRUE"="#0570b0")) +
  facet_nested(rows = vars(LongRead_MAG_in_final_cluster),
               scales = "free_y",
               space = "free_y",
               switch = "both") +
  scale_x_discrete(expand = c(0, 0), position = "top") +
  theme_classic() +
  theme_RE +
  theme_tiles +
  xlab("S40")

fig3_shortread20 <- df %>%
  ggplot(aes(x = "", y = reorder(OTU, Abundance), 
             fill = ShortRead20_MAG_in_final_cluster)) +
  geom_tile(colour = "white") +
  scale_fill_manual(values=c("FALSE"="white","TRUE"="#74a9cf")) +
  facet_nested(rows = vars(LongRead_MAG_in_final_cluster),
               scales = "free_y",
               space = "free_y",
               switch = "both") +
  scale_x_discrete(expand = c(0, 0), position = "top") +
  theme_classic() +
  theme_RE +
  theme_tiles +
  xlab("S20")

fig3_hybrid <- df %>%
  ggplot(aes(x = "", y = reorder(OTU, Abundance), 
             fill = Hybrid_MAG_in_final_cluster)) +
  geom_tile(colour = "white") +
  scale_fill_manual(values=c("FALSE"="white","TRUE"="#fc8d59")) +
  facet_nested(rows = vars(LongRead_MAG_in_final_cluster),
               scales = "free_y",
               space = "free_y",
               switch = "both") +
  scale_x_discrete(expand = c(0, 0), position = "top") +
  theme_classic() +
  theme_RE +
  theme_tiles +
  xlab("Hy")


fig3_longread <- df %>%
  ggplot(aes(x = "", y = reorder(OTU, Abundance), 
             fill = LongRead_MAG_in_final_cluster)) +
  geom_tile(colour = "white") +
  scale_fill_manual(values=c("MAGs not found by long-reads"="white",
                             "MAGs in all strategies"="#969696")) +
  facet_nested(rows = vars(LongRead_MAG_in_final_cluster),
               scales = "free_y",
               space = "free_y",
               switch = "both") +
  scale_x_discrete(expand = c(0, 0), position = "top") +
  theme_classic() +
  theme_RE +
  theme_tiles +
  xlab("LR")


fig3_completeness <- df %>%
  group_by(OTU) %>%
  summarise(completeness = unique(completeness),
            LongRead_MAG_in_final_cluster = unique(LongRead_MAG_in_final_cluster),
            mean_rel_abund = mean(Abundance)) %>%
  ungroup() %>%
  ggplot(aes(x = completeness, y = reorder(OTU, mean_rel_abund))) +
  geom_col(position = "stack") +
  facet_nested(rows = vars(LongRead_MAG_in_final_cluster),
               scales = "free_y",
               space = "free_y",
               switch = "both") +
  scale_x_continuous(expand = c(0, 0), 
                     position = "top", 
                     breaks = c(70, 80, 90)) +
  coord_cartesian(xlim = c(70, 100)) + 
  geom_vline(xintercept = c(80, 90, 100), colour = "darkgrey", alpha = 0.5) +
  theme_classic() +
  theme_RE +
  theme(
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_blank(),
    axis.title.x = element_markdown(size = 12),
    axis.title.y = element_blank(),
    axis.line.x = element_blank(),
  ) +
  xlab("Comp")


fig3_contamination <- df %>%
  group_by(OTU) %>%
  summarise(contamination = unique(contamination),
            LongRead_MAG_in_final_cluster = unique(LongRead_MAG_in_final_cluster),
            mean_rel_abund = mean(Abundance)) %>%
  ungroup() %>%
  ggplot(aes(x = contamination, y = reorder(OTU, mean_rel_abund))) +
  geom_col(position = "stack") +
  facet_nested(rows = vars(LongRead_MAG_in_final_cluster),
               scales = "free_y",
               space = "free_y",
               switch = "both") +
  scale_x_continuous(expand = c(0, 0), 
                     position = "top", 
                     breaks = c(0, 3, 6, 9)) +
  coord_cartesian(xlim = c(0, 10)) + 
  geom_vline(xintercept = c(3, 6, 9), colour = "darkgrey", alpha = 0.5) +
  theme_classic() +
  theme_RE +
  theme(
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_blank(),
    axis.title.x = element_markdown(size = 12),
    axis.title.y = element_blank(),
    axis.line.x = element_blank(),
  ) +
  xlab("Cont")


fig3_mean_rel_abun <- df %>%
  group_by(OTU) %>%
  summarise(LongRead_MAG_in_final_cluster = unique(LongRead_MAG_in_final_cluster),
            mean_rel_abund = mean(Abundance)) %>%
  ungroup() %>%
  ggplot(aes(x = mean_rel_abund * 100, y = reorder(OTU, mean_rel_abund))) +
  geom_col(position = "stack") +
  facet_nested(rows = vars(LongRead_MAG_in_final_cluster),
               scales = "free_y",
               space = "free_y",
               switch = "both") +
  scale_x_continuous(expand = c(0, 0), 
                     position = "top",
                     limits = c(0, 13),
                     breaks = c(0, 5, 10)) +
  theme_classic() +
  theme_RE +
  theme(
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_blank(),
    axis.title.x = element_markdown(size = 12),
    axis.title.y = element_blank(),
    axis.line.x = element_blank(),
  ) +
  xlab("M-RA")


design <- "111123333333333333456789a"

guide_area() + fig3_taxa_labels + fig3a + fig3_mean_rel_abun + fig3_shortread40 + fig3_shortread20 + fig3_hybrid + fig3_longread + fig3_completeness + fig3_contamination +
  plot_layout(design = design,
              heights = c(1, 6),
              guides = "collect") &
  theme(plot.margin = margin(0, 0, 0, 0),
        legend.margin = margin(0, -10, 0, 0),
  )


ggsave("figures/Figure_4.png", width = 12, height = 16, units = "in")
