library(tidyverse)
library(ggpubr)
library(viridis)
library(cowplot)
library(phyloseq)
library(patchwork)

phyloseq_file = 'data/18S/phyloseqs.RData' 

# Load phyloseq objects
load(phyloseq_file)

# Transform to relative abundance
phyloseq_mb = transform_sample_counts(phyloseq_mb, function(x)100* x / sum(x))
phyloseq_abutab = transform_sample_counts(phyloseq_abutab, function(x)100* x / sum(x))

# Remove unknown genera
phyloseq_abutab = subset_taxa(phyloseq_abutab, !taxon_genus=="")
phyloseq_mb = subset_taxa(phyloseq_mb, !Genus=="")

# Aggregate on genus level
phyloseq_seqtab_glom = tax_glom(phyloseq_mb, taxrank="Genus", NArm = FALSE)
phyloseq_abutab_glom = tax_glom(phyloseq_abutab, taxrank="taxon_genus",NArm = FALSE)

# Find top 15 genus from metabarcoding

top15A_counts<- sort(taxa_sums(phyloseq_seqtab_glom), decreasing = TRUE)[1:15]
top15_A <- names(sort(taxa_sums(phyloseq_seqtab_glom), decreasing = TRUE))[1:15]
top15_genus_mb = tax_table(phyloseq_seqtab_glom)[top15_A,]

# Find top 15 genus from microscopy

top15A_counts<- sort(taxa_sums(phyloseq_abutab_glom), decreasing = TRUE)[1:15]
top15_A <- names(sort(taxa_sums(phyloseq_abutab_glom), decreasing = TRUE))[1:15]
top15_genus_mc = tax_table(phyloseq_abutab_glom)[top15_A,]

# Bind top taxa tables together
top_taxa = rbind(top15_genus_mb[,8], top15_genus_mc[,7])

genus = c(unique(top_taxa[,1]))

taxa_in_common = top15_genus_mb[,8][top15_genus_mb[,8] %in% top15_genus_mc[,7]]


# Melt tables and calculate genus-average

psmelt_abutab = psmelt(phyloseq_abutab_glom) %>%
  filter(taxon_genus %in% genus) %>%
  filter(!Abundance == 0) %>%
  dplyr::group_by(taxon_genus) %>%
  dplyr::summarise("mean" = mean(Abundance),
                   "sd" = sd(Abundance),
                   "n" = n()) %>%
  mutate(method = "Utermöhl") %>%
  dplyr::rename(Genus = taxon_genus)
  
psmelt_seqtab = psmelt(phyloseq_seqtab_glom) %>%
  filter(Genus %in% genus) %>%
  filter(!Abundance == 0) %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarise("mean" = mean(Abundance),
                   "sd" = sd(Abundance),
                   "n" = n()) %>%
  mutate(method = "Metabarcoding")

# Bind tables
psmelt = rbind(psmelt_abutab, psmelt_seqtab)

# Plots

relative_abundance.p = ggplot(psmelt, aes(x = fct_reorder(Genus, desc(as.numeric(mean))), y = mean, fill = method)) +
  geom_bar(position="dodge", 
           stat="identity", 
           color = "black",
           linewidth = 0.25) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=0,
                position=position_dodge(.9), linewidth = 0.25) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none",
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  ylab("Relative abundance (%)") +
  xlab("") +
  scale_fill_manual(values=c('white', '#999999')) 

frequency.p = ggplot(psmelt, aes(x = fct_reorder(Genus, desc(as.numeric(mean))), y = n, fill = method)) +
  geom_bar(position="dodge", 
           stat="identity", 
           color = "black",
           linewidth = 0.25) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="top",
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  ylab("Frequency of occurrence (n)") +
  xlab("") +
  scale_fill_manual(values=c('white', '#999999')) 

# Extract the legend from one of the plots
legend <- get_legend(
  frequency.p + guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)


# Plots in grid

plot = plot_grid(relative_abundance.p, 
          frequency.p + theme(legend.position = "none"), 
          legend, 
          labels = c("a)", "b)"), 
          label_x = 0.93,
          label_y = 0.97,
          ncol = 1,
          rel_heights = c(1, 1, .1))

# Save plot

ggsave("fig4_top_genera.png", 
       width = 7,
       height = 7,
       plot = plot, 
       device = "png", 
       path = "plots/microscopy_vs_barcoding/",
       bg = "white")


### Prepare heatmap

# Melt tables and identify samples with presence
psmelt_mb_heatmap = psmelt(phyloseq_seqtab_glom) %>%
  filter(Genus %in% genus) %>%
  dplyr::select(NGI_sample_ID_18S,
                sampling_date,
                station_name_alt,
                sea_basin,
                day_of_year,
                season,
                Genus,
                Abundance) %>%
  dplyr::rename(Abundance_mb = Abundance) %>%
  mutate(Abundance_mb = ifelse(is.na(Abundance_mb), 0, Abundance_mb)) %>%
  mutate("in_mb" = ifelse(Abundance_mb > 0, "Yes", "No"))

psmelt_um_heatmap = psmelt(phyloseq_abutab_glom) %>%
  filter(taxon_genus %in% genus) %>%
  dplyr::select(NGI_sample_ID_18S,
                sampling_date,
                station_name_alt,
                sea_basin,
                day_of_year,
                season,
                taxon_genus,
                Abundance) %>%
  dplyr::rename(Genus = taxon_genus) %>%
  dplyr::rename(Abundance_um = Abundance) %>%
  mutate(Abundance_um = ifelse(is.na(Abundance_um), 0, Abundance_um)) %>%
  mutate("in_um" = ifelse(Abundance_um > 0, "Yes", "No"))

# Join table and wrangle data
psmelt_heatmap = psmelt_mb_heatmap %>%
  full_join(psmelt_um_heatmap) %>%
  mutate(Abundance_um = ifelse(is.na(Abundance_um), 0, Abundance_um)) %>%
  mutate(Abundance_mb = ifelse(is.na(Abundance_mb), 0, Abundance_mb)) %>%
  mutate("in_um" = ifelse(Abundance_um > 0, "Yes", "No")) %>%
  mutate("in_mb" = ifelse(Abundance_mb > 0, "Yes", "No")) %>%
  mutate("Presence" = ifelse(in_mb == "Yes" & in_um == "Yes", "MB & UM", NA)) %>%
  mutate("Presence_mb" = ifelse(in_mb == "Yes" & in_um != "Yes", "MB only", NA))%>%
  mutate("Presence_um" = ifelse(in_mb != "Yes" & in_um == "Yes", "UM only", NA)) %>%
  mutate(Presence = coalesce(Presence, Presence_mb, Presence_um)) %>%
  dplyr::select(-Presence_mb, 
                -Presence_um,
                -Abundance_mb,
                -Abundance_um,
                -in_mb,
                -in_um) %>%
  mutate(sample_name = paste(sampling_date, station_name_alt)) %>%
  mutate(sample_name = factor(sample_name)) %>%
  mutate(basin_name = sea_basin) %>%
  mutate(basin_name = dplyr::recode(basin_name, 'Bothnian Bay' = 1, 'Bothnian Sea' = 2, 'Baltic Proper' = 3, 'Kattegat' = 4, 'Skagerrak' = 5)) %>%
  mutate(sample_name = fct_reorder(sample_name, as.Date(sampling_date))) %>%
  mutate(Genus = factor(Genus))

psmelt_heatmap %>%
  group_by(Presence) %>%
  tally()

# Sum of samples
sum_mb = psmelt_heatmap %>%
  dplyr::group_by(Genus) %>%
  filter(Presence == "MB & UM" | Presence == "MB only") %>%
  dplyr::tally() %>%
  mutate("Method" = "MB (\U1D62F)")

missing_mb = genus[!genus %in% sum_mb$Genus]

sum_mb = sum_mb %>%
  add_row(Genus = missing_mb, n = 0, Method = "MB (\U1D62F)")

sum_um = psmelt_heatmap %>%
  dplyr::group_by(Genus) %>%
  filter(Presence == "MB & UM" | Presence == "UM only") %>%
  dplyr::tally() %>%
  mutate("Method" = "UM (\U1D62F)")

missing_um = genus[!genus %in% sum_um$Genus]

sum_um = sum_um %>%
  add_row(Genus = missing_um, n = 0, Method = "UM (\U1D62F)")

sum = rbind(sum_mb, sum_um)

tot = sum %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarise(tot = sum(n))

psmelt_heatmap = psmelt_heatmap %>%
  left_join(tot) %>%
  mutate(Genus = fct_reorder(Genus, tot, .desc = TRUE))
  
sum = sum %>%
  left_join(tot) %>%
  mutate(Genus = fct_reorder(Genus, tot, .desc = TRUE)) %>%
  mutate(label = "Sum")

psmelt_heatmap %>%
  filter(!is.na(Presence)) %>%
  group_by(Presence) %>%
  summarise(n = n()) %>%
  mutate(percent = n/sum(n)*100)

# Plot heatmap
heatmap <- ggplot(psmelt_heatmap, aes(Genus, sample_name, label = season)) +
  geom_tile(aes(fill = Presence)) +
  geom_vline(xintercept = seq(0.5, nlevels(psmelt_heatmap$Genus) + 0.5), color = "black", linetype = "solid", linewidth = .1) +  # Add vertical lines
  theme(axis.text.x = element_text(margin = unit(c(0, 0, 0, 0), "mm"), angle = 45, vjust = 1, hjust = 1),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "grey25", fill = NA, linewidth = .5),
        strip.background = element_rect(fill = "white", color = "grey25", linewidth = .5),
        panel.grid.major = element_blank(), 
        panel.background = element_blank(),
        panel.spacing.x = unit(2, "lines"),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        axis.text.y=element_blank(),
        text = element_text(size = 9),
        panel.spacing = unit(0.05,'lines'),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.margin=margin(t = 0, unit='cm'),
        legend.box.margin=margin(-10,-10,0,0),
        legend.key.height = unit(.1, 'cm'),
        legend.text = element_text(angle = 0, hjust = 1, vjust = 0.5)
  ) +
  facet_grid(rows = vars(sea_basin), scales = "free", space = "free") +
  ylab("") +
  xlab("") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  # scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3"))
  scale_fill_viridis(discrete = TRUE, option = "D", alpha = 0.7, na.value = "grey90")

# Table

table = ggplot(sum, aes(x = Genus, y = Method)) +
  geom_tile(fill = NA, color = "black") +
  geom_text(aes(x = Genus, y = Method, label = n), color = "grey30", size = 2.5) +
  labs(y = "", x = NULL) +
  facet_grid(rows=vars(label)) +
  theme_minimal() +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  theme(axis.line = element_blank(), 
        plot.margin = margin(0, 0, 0, 0, "pt"),
        # panel.grid.major.x = element_line(color = "red",
        #                                   linewidth = 0.5,
        #                                   linetype = 2),
        text = element_text(size = 9),
        panel.border = element_rect(color = "grey25", fill = NA, linewidth = .5),
        axis.text.y=element_text(size = 7, margin = margin(r = 0)),
        axis.ticks = element_blank(), 
        strip.background = element_rect(fill = "white", color = "grey25", linewidth = .5),
        axis.text.x = element_blank(),
        panel.spacing = unit(0, "mm"),
        panel.grid = element_blank())

heatmap_table = table / heatmap + plot_layout(heights = c(1, 20))

# Save plot
ggsave("fig4_heatmap.png", 
       width = 5,
       height = 7,
       dpi = 300,
       plot = heatmap_table, 
       device = "png", 
       path = "plots/microscopy_vs_barcoding/",
       bg = "white")

# Save PDF
ggsave(
  plot = heatmap_table,
  path = "plots/microscopy_vs_barcoding/pdf",
  filename = "fig4_heatmap.pdf",
  device = "pdf",
  width = 5,
  height = 7,
  bg = "white"
)
