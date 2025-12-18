library(tidyverse)
library(phyloseq)
library(cowplot)

phyloseq_file = 'data/18S/phyloseqs.RData' 

# Load phyloseq objects
load(phyloseq_file)

# Aggregate on class level
phyloseq_seqtab_glom = tax_glom(phyloseq_mb, taxrank="Class", NArm = FALSE)
phyloseq_abutab_glom = tax_glom(phyloseq_abutab, taxrank="Class_MB", NArm = FALSE)
phyloseq_voltab_glom = tax_glom(phyloseq_voltab, taxrank="Class_MB", NArm = FALSE)
phyloseq_carbtab_glom = tax_glom(phyloseq_carbtab, taxrank="Class_MB", NArm = FALSE)

# Transform to relative abundance
phyloseq_seqtab_glom = transform_sample_counts(phyloseq_seqtab_glom, function(x)100* x / sum(x))
phyloseq_abutab_glom = transform_sample_counts(phyloseq_abutab_glom, function(x)100* x / sum(x))
phyloseq_voltab_glom = transform_sample_counts(phyloseq_voltab_glom, function(x)100* x / sum(x))
phyloseq_carbtab_glom = transform_sample_counts(phyloseq_carbtab_glom, function(x)100* x / sum(x))

# Melt tables and calculate average abundance per sea basin

# Melt mc abundance
psmelt_abutab = psmelt(phyloseq_abutab_glom) %>%
  # filter(!Abundance == 0) %>%
  dplyr::group_by(Class_MB, sea_basin) %>%
  dplyr::summarise("Abundance" = mean(Abundance)) %>%
  mutate(method = "Utermöhl:\nabundance") %>%
  dplyr::rename(Class = Class_MB) %>%
  mutate(basin_name = sea_basin) %>%
  mutate(basin_name = dplyr::recode(basin_name, 'Bothnian Bay' = 'I', 'Bothnian Sea' = 'II', 'Baltic Proper' = 'III', 'Kattegat' = 'IV', 'Skagerrak' = 'V'))

# Melt mc biovolume concentration
psmelt_voltab = psmelt(phyloseq_voltab_glom) %>%
  # filter(!Abundance == 0) %>%
  dplyr::group_by(Class_MB, sea_basin) %>%
  dplyr::summarise("Abundance" = mean(Abundance)) %>%
  mutate(method = "Utermöhl:\nbiovolume concentration") %>%
  dplyr::rename(Class = Class_MB) %>%
  mutate(basin_name = sea_basin) %>%
  mutate(basin_name = dplyr::recode(basin_name, 'Bothnian Bay' = 'I', 'Bothnian Sea' = 'II', 'Baltic Proper' = 'III', 'Kattegat' = 'IV', 'Skagerrak' = 'V'))

# Melt mc carbon concentration
psmelt_carbtab = psmelt(phyloseq_carbtab_glom) %>%
  # filter(!Abundance == 0) %>%
  dplyr::group_by(Class_MB, sea_basin) %>%
  dplyr::summarise("Abundance" = mean(Abundance)) %>%
  mutate(method = "Utermöhl:\ncarbon concentration") %>%
  dplyr::rename(Class = Class_MB) %>%
  mutate(basin_name = sea_basin) %>%
  mutate(basin_name = dplyr::recode(basin_name, 'Bothnian Bay' = 'I', 'Bothnian Sea' = 'II', 'Baltic Proper' = 'III', 'Kattegat' = 'IV', 'Skagerrak' = 'V'))

# Melt mb reads
psmelt_seqtab = psmelt(phyloseq_seqtab_glom) %>%
  # filter(!Abundance == 0) %>%
  dplyr::group_by(Class, sea_basin) %>%
  dplyr::summarise("Abundance" = mean(Abundance)) %>%
  mutate(method = "Metabarcoding:\nread counts") %>%
  mutate(basin_name = sea_basin) %>%
  mutate(basin_name = dplyr::recode(basin_name, 'Bothnian Bay' = 'I', 'Bothnian Sea' = 'II', 'Baltic Proper' = 'III', 'Kattegat' = 'IV', 'Skagerrak' = 'V'))

# Check sums = 100
psmelt_seqtab %>%
  ungroup %>%
  dplyr::group_by(sea_basin) %>%
  dplyr::summarise("sum" = sum(Abundance))

# Bind tables and format for plot
psmelt = rbind(psmelt_abutab, 
               psmelt_seqtab,
               psmelt_voltab,
               psmelt_carbtab) %>%
  mutate(method = factor(method, levels = c("Metabarcoding:\nread counts", "Utermöhl:\nabundance", "Utermöhl:\ncarbon concentration", "Utermöhl:\nbiovolume concentration"))) %>%
  mutate(Class = gsub("ophyceae", "", Class)) %>%
  mutate(basin_name = dplyr::recode(basin_name, 'Bothnian Bay' = 'I', 'Bothnian Sea' = 'II', 'Baltic Proper' = 'III', 'Kattegat' = 'IV', 'Skagerrak' = 'V'))

# Summary for all data

# Metabarcoding reads
psmelt(phyloseq_seqtab_glom) %>%
  dplyr::group_by(Class) %>%
  dplyr::summarise("Abundance" = mean(Abundance)) %>%
  arrange(desc(Abundance))

# abundance
psmelt(phyloseq_abutab_glom) %>%
  dplyr::group_by(Class_MB) %>%
  dplyr::summarise("Abundance" = mean(Abundance)) %>%
  arrange(desc(Abundance))

# Carbon concentration
psmelt(phyloseq_carbtab_glom) %>%
  dplyr::group_by(Class_MB) %>%
  dplyr::summarise("Abundance" = mean(Abundance)) %>%
  arrange(desc(Abundance))

# Biovolume concentration
psmelt(phyloseq_voltab_glom) %>%
  dplyr::group_by(Class_MB) %>%
  dplyr::summarise("Abundance" = mean(Abundance))%>%
  arrange(desc(Abundance))

# Color palettes
c23 <- c(
  "dodgerblue2", 
  "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "skyblue2", 
  "#FB9A99", # lt pink
  "#CAB2D6", # lt purple
  "palegreen2",
  "#FDBF6F", # lt orange
  "gray70", 
  "khaki2",
  "maroon", 
  "orchid1", 
  "deeppink1", 
  "blue1", 
  "steelblue4",
  "darkturquoise", 
  "green1", 
  "yellow4", 
  "yellow3",
  "darkorange4", 
  "brown"
)

# Plot
relative_abundance.p = ggplot(psmelt, aes(x = fct_reorder(basin_name, as.numeric(sea_basin)), y = Abundance, fill = Class)) +
  geom_bar(position="stack", 
           stat="identity", 
           color = "black",
           linewidth = 0.25) +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.key.height = unit(.2, 'cm'),
        legend.spacing.y = unit(.1, 'cm'),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank()) +
  ylab("Relative proportion (%)") +
  # xlab("") +
  guides(fill = guide_legend(byrow = TRUE)) +
  facet_grid(cols = vars(method)) +
  scale_fill_manual(values = alpha(c23, 0.6)) +
  xlab("Sea basin")

# Save plot
ggsave("fig6_classes_abundance.png",
       width = 8,
       height = 5,
       plot = relative_abundance.p,
       dpi = 300,
       device = "png",
       path = "plots/microscopy_vs_barcoding/",
       bg = "white")

# Save PDF
ggsave(
  plot = relative_abundance.p,
  path = "plots/microscopy_vs_barcoding/pdf",
  filename = "fig6_classes.pdf",
  device = "pdf",
  width = 8,
  height = 5,
  bg = "white"
)

# Save table
psmelt_table <- psmelt %>%
  mutate(method = gsub("\n", "", method))

write_tsv(psmelt_table, "Fig6.txt")

# Seasonal plots for the high frequency stations

# Melt tables and select stations for seasonal plot
psmelt_abutab_season = psmelt(phyloseq_abutab_glom) %>%
  filter(station_name == "SL\x80GG\x85" | station_name == "B7" | station_name == "B1") %>%
  mutate(station_name = gsub("SL\x80GG\x85", "Släggö", station_name)) %>%
  mutate(station_name = gsub("B7", "B3 / B7", station_name)) %>%
  # dplyr::group_by(Class_MB, sea_basin) %>%
  # dplyr::summarise("Abundance" = mean(Abundance)) %>%
  mutate(method = "Utermöhl:\nabundance") %>%
  dplyr::rename(Class = Class_MB) %>%
  select(-station_date, -Flagellates, -Unicell)

psmelt_voltab_season = psmelt(phyloseq_voltab_glom) %>%
  filter(station_name == "SL\x80GG\x85" | station_name == "B7" | station_name == "B1") %>%
  mutate(station_name = gsub("SL\x80GG\x85", "Släggö", station_name)) %>%
  mutate(station_name = gsub("B7", "B3 / B7", station_name)) %>%
  # dplyr::group_by(Class_MB, sea_basin) %>%
  # dplyr::summarise("Abundance" = mean(Abundance)) %>%
  mutate(method = "Utermöhl:\nbiovolume concentration") %>%
  dplyr::rename(Class = Class_MB) %>%
  select(-station_date, -Flagellates, -Unicell)

psmelt_carbtab_season = psmelt(phyloseq_carbtab_glom) %>%
  filter(station_name == "SL\x80GG\x85" | station_name == "B7" | station_name == "B1") %>%
  mutate(station_name = gsub("SL\x80GG\x85", "Släggö", station_name)) %>%
  mutate(station_name = gsub("B7", "B3 / B7", station_name)) %>%
  # dplyr::group_by(Class_MB, sea_basin) %>%
  # dplyr::summarise("Abundance" = mean(Abundance)) %>%
  mutate(method = "Utermöhl:\ncarbon concentration") %>%
  dplyr::rename(Class = Class_MB) %>%
  select(-station_date, -Flagellates, -Unicell)

psmelt_seqtab_season = psmelt(phyloseq_seqtab_glom) %>%
  filter(station_name == "SL\x80GG\x85" | station_name == "B7" | station_name == "B1") %>%
  mutate(station_name = gsub("SL\x80GG\x85", "Släggö", station_name)) %>%
  mutate(station_name = gsub("B7", "B3 / B7", station_name)) %>%
  mutate(method = "Metabarcoding:\nread counts") %>%
  select(names(psmelt_abutab_season))

# Bind tables
psmelt_season = rbind(
  psmelt_carbtab_season,
  # psmelt_voltab_season,
  # psmelt_abutab_season,
  psmelt_seqtab_season
  ) %>%
  mutate(method = factor(method, levels = c("Metabarcoding:\nread counts", "Utermöhl:\nabundance", "Utermöhl:\ncarbon concentration", "Utermöhl:\nbiovolume concentration"))) %>%
  mutate(Class = gsub("ophyceae", "", Class))

# Plot
season.p = ggplot(psmelt_season, aes(x = sampling_date, y = Abundance, fill = Class)) +
  geom_bar(position="stack", 
           stat="identity", 
           color = "black",
           linewidth = 0.25) +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.key.height = unit(.2, 'cm'),
        legend.spacing.y = unit(.1, 'cm'),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank()) +
  ylab("Relative proportion (%)") +
  scale_x_discrete(labels=function(x) month.abb[lubridate::month(x)]) +
  guides(fill = guide_legend(byrow = TRUE)) +
  facet_grid(cols = vars(station_name), rows = vars(method), scales = "free_x") +
  scale_fill_manual(values = alpha(c23, 0.6)) +
  xlab("")

# Save plot
ggsave("figX_classes_season.png", 
       width = 10,
       height = 5,
       plot = season.p, 
       device = "png", 
       path = "plots/microscopy_vs_barcoding/",
       bg = "white")
