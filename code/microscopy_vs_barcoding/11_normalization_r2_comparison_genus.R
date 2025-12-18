library(tidyverse)
library(phyloseq)
library(cowplot)
library(patchwork)
library(viridis)
library(RColorBrewer)

# File paths
phyloseq_file = 'data/18S/phyloseqs.RData'

# Load phyloseq objects
load(phyloseq_file)

# Subset samples to identical samples
phyloseq_abutab_no_glom = prune_samples(sample_data(phyloseq_abutab)$NGI_sample_ID_18S %in%
                                          sample_data(phyloseq_spike)$NGI_sample_ID_18S, phyloseq_abutab)

phyloseq_carbtab_no_glom = prune_samples(sample_data(phyloseq_carbtab)$NGI_sample_ID_18S %in%
                                           sample_data(phyloseq_spike)$NGI_sample_ID_18S, phyloseq_carbtab)

phyloseq_voltab_no_glom = prune_samples(sample_data(phyloseq_voltab)$NGI_sample_ID_18S %in%
                                          sample_data(phyloseq_spike)$NGI_sample_ID_18S, phyloseq_voltab)

phyloseq_mb_no_glom = prune_samples(sample_data(phyloseq_mb)$NGI_sample_ID_18S %in%
                                      sample_data(phyloseq_spike)$NGI_sample_ID_18S, phyloseq_mb)

phyloseq_dna_no_glom = prune_samples(sample_data(phyloseq_dna)$NGI_sample_ID_18S %in%
                                       sample_data(phyloseq_spike)$NGI_sample_ID_18S, phyloseq_dna)

phyloseq_spike_no_glom = phyloseq_spike

# Aggregate TO genus level
phyloseq_mb = tax_glom(phyloseq_mb_no_glom, taxrank="Genus", NArm = FALSE)
phyloseq_spike = tax_glom(phyloseq_spike, taxrank="Genus", NArm = FALSE)
phyloseq_dna = tax_glom(phyloseq_dna_no_glom, taxrank="Genus", NArm = FALSE)
phyloseq_abutab = tax_glom(phyloseq_abutab_no_glom, taxrank="taxon_genus", NArm = FALSE)
phyloseq_carbtab = tax_glom(phyloseq_carbtab_no_glom, taxrank="taxon_genus", NArm = FALSE)
phyloseq_voltab = tax_glom(phyloseq_voltab_no_glom, taxrank="taxon_genus", NArm = FALSE)

# Transform to relative abundance
phyloseq_mb_rel = transform_sample_counts(phyloseq_mb, function(x)100* x / sum(x))
phyloseq_abutab_rel = transform_sample_counts(phyloseq_abutab, function(x)100* x / sum(x))
phyloseq_carbtab_rel = transform_sample_counts(phyloseq_carbtab, function(x)100* x / sum(x))
phyloseq_voltab_rel = transform_sample_counts(phyloseq_voltab, function(x)100* x / sum(x))

# Find top 50 genus from metabarcoding (relative abundance)
top50A_counts<- sort(taxa_sums(phyloseq_mb_rel), decreasing = TRUE)[1:50]
top50_A <- names(sort(taxa_sums(phyloseq_mb_rel), decreasing = TRUE))[1:50]
top50_genus_mb = phyloseq::tax_table(phyloseq_mb_rel)[top50_A,]

# Find top 50 genus from microscopy
top50A_counts<- sort(taxa_sums(phyloseq_abutab_rel), decreasing = TRUE)[1:50]
top50_B <- names(sort(taxa_sums(phyloseq_abutab_rel), decreasing = TRUE))[1:50]
top50_genus_mc = phyloseq::tax_table(phyloseq_abutab_rel)[top50_B,]

# Find top 20 genus from metabarcoding (relative abundance)
top15A_counts<- sort(taxa_sums(phyloseq_mb_rel), decreasing = TRUE)[1:20]
top15_A <- names(sort(taxa_sums(phyloseq_mb_rel), decreasing = TRUE))[1:20]
top15_genus_mb = phyloseq::tax_table(phyloseq_mb_rel)[top15_A,]

# Find top 20 genus from microscopy
top15A_counts<- sort(taxa_sums(phyloseq_abutab_rel), decreasing = TRUE)[1:20]
top15_B <- names(sort(taxa_sums(phyloseq_abutab_rel), decreasing = TRUE))[1:20]
top15_genus_mc = phyloseq::tax_table(phyloseq_abutab_rel)[top15_B,]


# Bind the tables
top_taxa = rbind(top50_genus_mb[,8], top50_genus_mc[,7])

top_15_taxa = rbind(top15_genus_mb[,8], top15_genus_mc[,7])


genus = c(unique(top_15_taxa[,1]))
genus <- genus[!is.na(genus)]
genus <- genus[!genus == ""]

taxa_in_common = top50_genus_mb[,8][top50_genus_mb[,8] %in% top50_genus_mc[,7]]

taxa_in_common <- genus

# Select only relevant genera
phyloseq_mb = subset_taxa(phyloseq_mb, Genus %in% taxa_in_common)
phyloseq_spike = subset_taxa(phyloseq_spike, Genus %in% taxa_in_common)
phyloseq_dna = subset_taxa(phyloseq_dna, Genus %in% taxa_in_common)
phyloseq_abutab = subset_taxa(phyloseq_abutab, taxon_genus %in% taxa_in_common)
phyloseq_carbtab = subset_taxa(phyloseq_carbtab, taxon_genus %in% taxa_in_common)
phyloseq_voltab = subset_taxa(phyloseq_voltab, taxon_genus %in% taxa_in_common)
phyloseq_mb_rel = subset_taxa(phyloseq_mb_rel, Genus %in% taxa_in_common)
phyloseq_abutab_rel = subset_taxa(phyloseq_abutab_rel, taxon_genus %in% taxa_in_common)
phyloseq_carbtab_rel = subset_taxa(phyloseq_carbtab_rel, taxon_genus %in% taxa_in_common)
phyloseq_voltab_rel = subset_taxa(phyloseq_voltab_rel, taxon_genus %in% taxa_in_common)

# Melt tables
psmelt_seqtab = psmelt(phyloseq_mb) %>%
  filter(!Abundance == 0) %>%
  mutate(method = "Metabarcoding") %>%
  mutate(abundance_type = "Reads") %>%
  select(NGI_sample_ID_18S,
         sea_basin,
         Genus,
         Abundance) %>%
  rename("mb_reads" = Abundance)

psmelt_seqtab_rel = psmelt(phyloseq_mb_rel) %>%
  filter(!Abundance == 0) %>%
  mutate(method = "Metabarcoding") %>%
  mutate(abundance_type = "Relative abundance") %>%
  select(NGI_sample_ID_18S,
         sea_basin,
         Genus,
         Abundance) %>%
  rename("mb_relative" = Abundance)

psmelt_spiketab = psmelt(phyloseq_spike) %>%
  filter(!Abundance == 0) %>%
  mutate(method = "Metabarcoding") %>%
  mutate(abundance_type = "Spike-normalized") %>%
  select(NGI_sample_ID_18S,
         sea_basin,
         Genus,
         Abundance) %>%
  rename("mb_spike" = Abundance)

psmelt_dnatab = psmelt(phyloseq_dna) %>%
  filter(!Abundance == 0) %>%
  mutate(method = "Metabarcoding") %>%
  mutate(abundance_type = "DNA-normalized") %>%
  select(NGI_sample_ID_18S,
         sea_basin,
         Genus,
         Abundance) %>%
  rename("mb_dna" = Abundance)

psmelt_abutab = psmelt(phyloseq_abutab) %>%
  filter(!Abundance == 0) %>%
  mutate(method = "Microscopy") %>%
  mutate(abundance_type = "Cell counts") %>%
  dplyr::rename(Genus = taxon_genus) %>%
  select(NGI_sample_ID_18S,
         sea_basin,
         Genus,
         Abundance) %>%
  rename("mc_counts" = Abundance)

psmelt_carbtab = psmelt(phyloseq_carbtab) %>%
  filter(!Abundance == 0) %>%
  mutate(method = "Microscopy") %>%
  mutate(abundance_type = "Carbon concentration") %>%
  dplyr::rename(Genus = taxon_genus) %>%
  select(NGI_sample_ID_18S,
         sea_basin,
         Genus,
         Abundance) %>%
  rename("mc_carbon" = Abundance)

psmelt_voltab = psmelt(phyloseq_voltab) %>%
  filter(!Abundance == 0) %>%
  mutate(method = "Microscopy") %>%
  mutate(abundance_type = "Biovolume concentration") %>%
  dplyr::rename(Genus = taxon_genus) %>%
  select(NGI_sample_ID_18S,
         sea_basin,
         Genus,
         Abundance) %>%
  rename("mc_volume" = Abundance)

psmelt_abutab_rel = psmelt(phyloseq_abutab_rel) %>%
  filter(!Abundance == 0) %>%
  mutate(method = "Microscopy") %>%
  mutate(abundance_type = "Relative abundance") %>%
  dplyr::rename(Genus = taxon_genus) %>%
  select(NGI_sample_ID_18S,
         sea_basin,
         Genus,
         Abundance) %>%
  rename("mc_relative" = Abundance)

psmelt_carbtab_rel = psmelt(phyloseq_carbtab_rel) %>%
  filter(!Abundance == 0) %>%
  mutate(method = "Microscopy") %>%
  mutate(abundance_type = "Carbon relative abundance") %>%
  dplyr::rename(Genus = taxon_genus) %>%
  select(NGI_sample_ID_18S,
         sea_basin,
         Genus,
         Abundance) %>%
  rename("carbon_relative" = Abundance)

psmelt_voltab_rel = psmelt(phyloseq_voltab_rel) %>%
  filter(!Abundance == 0) %>%
  mutate(method = "Microscopy") %>%
  mutate(abundance_type = "Biovolume relative abundance") %>%
  dplyr::rename(Genus = taxon_genus) %>%
  select(NGI_sample_ID_18S,
         sea_basin,
         Genus,
         Abundance) %>%
  rename("biovol_relative" = Abundance)

# Check sums = 100
psmelt_abutab_rel %>%
  ungroup %>%
  dplyr::group_by(NGI_sample_ID_18S) %>%
  dplyr::summarise("sum" = sum(mc_relative)) %>%
  dplyr::summarise("max" = max(sum),
                   "min" = min(sum))

psmelt_seqtab_rel %>%
  ungroup %>%
  dplyr::group_by(NGI_sample_ID_18S) %>%
  dplyr::summarise("sum" = sum(mb_relative)) %>%
  dplyr::summarise("max" = max(sum),
                   "min" = min(sum))

psmelt_carbtab_rel %>%
  ungroup %>%
  dplyr::group_by(NGI_sample_ID_18S) %>%
  dplyr::summarise("sum" = sum(carbon_relative)) %>%
  dplyr::summarise("max" = max(sum),
                   "min" = min(sum))

psmelt_voltab_rel %>%
  ungroup %>%
  dplyr::group_by(NGI_sample_ID_18S) %>%
  dplyr::summarise("sum" = sum(biovol_relative)) %>%
  dplyr::summarise("max" = max(sum),
                   "min" = min(sum))

# Join melted tables
psmelt_all = psmelt_seqtab_rel %>%
  left_join(psmelt_abutab_rel) %>%
  left_join(psmelt_abutab) %>%
  left_join(psmelt_carbtab) %>%
  left_join(psmelt_voltab) %>%
  left_join(psmelt_seqtab) %>%
  left_join(psmelt_spiketab) %>%
  left_join(psmelt_dnatab) %>%
  left_join(psmelt_carbtab_rel) %>%
  left_join(psmelt_voltab_rel) %>%
  filter(!is.na(mc_counts))

# Add log-transformed versions (log1p to handle zeros)
psmelt_all <- psmelt_all %>%
  mutate(
    mb_reads_log     = log1p(mb_reads),
    mb_spike_log     = log1p(mb_spike),
    mb_dna_log       = log1p(mb_dna),
    mb_relative_log  = log1p(mb_relative),
    mc_counts_log    = log1p(mc_counts),
    mc_carbon_log    = log1p(mc_carbon),
    mc_volume_log    = log1p(mc_volume),
    mc_relative_log  = log1p(mc_relative),
    carbon_relative_log = log1p(carbon_relative),
    biovol_relative_log = log1p(biovol_relative)
  )

# Define comparisons
comparisons <- list(
  mb_reads  = c("mc_counts", "mc_carbon", "mc_volume"),
  mb_spike  = c("mc_counts", "mc_carbon", "mc_volume"),
  mb_dna    = c("mc_counts", "mc_carbon", "mc_volume"),
  mb_relative = c("mc_relative", "carbon_relative", "biovol_relative")
)

# Define log-transformed comparisons
comparisons_log <- list(
  mb_reads_log     = c("mc_counts_log", "mc_carbon_log", "mc_volume_log"),
  mb_spike_log     = c("mc_counts_log", "mc_carbon_log", "mc_volume_log"),
  mb_dna_log       = c("mc_counts_log", "mc_carbon_log", "mc_volume_log"),
  mb_relative_log  = c("mc_relative_log", "carbon_relative_log", "biovol_relative_log")
)

# ---- Helper functions ----
compute_correlations <- function(df, comparisons) {
  map_dfr(names(comparisons), function(mb_var) {
    map_dfr(comparisons[[mb_var]], function(mc_var) {
      cor_val <- suppressWarnings(cor(df[[mb_var]], df[[mc_var]],
                                      use = "complete.obs", method = "spearman"))
      tibble(mb_var = mb_var, mc_var = mc_var, correlation = cor_val)
    })
  })
}

compute_r2 <- function(df, comparisons) {
  map_dfr(names(comparisons), function(mb_var) {
    map_dfr(comparisons[[mb_var]], function(mc_var) {
      fit <- tryCatch(lm(df[[mb_var]] ~ df[[mc_var]]), error = function(e) NULL)
      tibble(
        mb_var = mb_var,
        mc_var = mc_var,
        r2 = if (!is.null(fit)) summary(fit)$r.squared else NA_real_
      )
    })
  })
}

# ---- Compute by genus ----
cor_genus <- psmelt_all %>%
  group_by(Genus) %>%
  group_modify(~ compute_correlations(.x, comparisons = comparisons_log)) %>%
  ungroup()

r2_genus <- psmelt_all %>%
  group_by(Genus) %>%
  group_modify(~ compute_r2(.x, comparisons = comparisons_log)) %>%
  ungroup()

# ---- Prepare plotting data ----
mc_mapping <- c(
  mc_counts = "Counts", mc_relative = "Counts",
  mc_carbon = "Carbon", carbon_relative = "Carbon",
  mc_volume = "Biovolume", biovol_relative = "Biovolume"
)

mb_labels <- c(
  mb_reads    = "Reads",
  mb_spike    = "Spike-norm",
  mb_dna      = "DNA-norm",
  mb_relative = "Relative (%)"
)

cor_genus_plot <- cor_genus %>%
  filter(
    (mb_var %in% c("mb_reads_log", "mb_spike_log", "mb_dna_log") & mc_var %in% c("mc_counts_log", "mc_carbon_log", "mc_volume_log")) |
      (mb_var == "mb_relative" & mc_var %in% c("mc_relative_log", "carbon_relative_log", "biovol_relative_log"))
  ) %>%
  mutate(
    mc_group = mc_mapping[mc_var],
    mb_label = mb_labels[mb_var],
    mc_group = factor(mc_group, levels = c("Counts", "Carbon", "Biovolume"))
  )

r2_genus_plot <- r2_genus %>%
  filter(
    (mb_var %in% c("mb_reads", "mb_spike", "mb_dna") & mc_var %in% c("mc_counts", "mc_carbon", "mc_volume")) |
      (mb_var == "mb_relative" & mc_var %in% c("mc_relative", "carbon_relative", "biovol_relative"))
  ) %>%
  mutate(
    mc_group = mc_mapping[mc_var],
    mb_label = mb_labels[mb_var],
    mc_group = factor(mc_group, levels = c("Counts", "Carbon", "Biovolume"))
  )

# Ensure facet order (for both correlation and R² plots)
panel_order <- c("Reads", "Spike-norm", "DNA-norm", "Relative (%)")

cor_genus_plot <- cor_genus_plot %>%
  mutate(
    mb_label = factor(mb_label, levels = panel_order),
    Genus = factor(Genus, levels = rev(sort(unique(Genus))))  # reverse for top-to-bottom A→Z
  )

r2_genus_plot <- r2_genus_plot %>%
  mutate(
    mb_label = factor(mb_label, levels = panel_order),
    Genus = factor(Genus, levels = rev(sort(unique(Genus))))  # reverse for top-to-bottom A→Z
  )

### LOG TRANSFORM ###

# ---- 1. Log-transform the relevant columns ----
psmelt_all <- psmelt_all %>%
  mutate(
    # Metabarcoding
    mb_reads_log     = log1p(mb_reads),
    mb_spike_log     = log1p(mb_spike),
    mb_dna_log       = log1p(mb_dna),
    mb_relative_log  = log1p(mb_relative),
    # Microscopy
    mc_counts_log    = log1p(mc_counts),
    mc_carbon_log    = log1p(mc_carbon),
    mc_volume_log    = log1p(mc_volume),
    mc_relative_log      = log1p(mc_relative),
    carbon_relative_log  = log1p(carbon_relative),
    biovol_relative_log  = log1p(biovol_relative)
  )

# ---- 2. Define comparisons (log-transformed) ----
comparisons_log <- list(
  mb_reads_log     = c("mc_counts_log", "mc_carbon_log", "mc_volume_log"),
  mb_spike_log     = c("mc_counts_log", "mc_carbon_log", "mc_volume_log"),
  mb_dna_log       = c("mc_counts_log", "mc_carbon_log", "mc_volume_log"),
  mb_relative_log  = c("mc_relative_log", "carbon_relative_log", "biovol_relative_log")
)

# ---- 3. Helper functions ----
compute_correlations <- function(df, comparisons) {
  map_dfr(names(comparisons), function(mb_var) {
    map_dfr(comparisons[[mb_var]], function(mc_var) {
      cor_val <- suppressWarnings(
        cor(df[[mb_var]], df[[mc_var]], use = "complete.obs", method = "spearman")
      )
      tibble(mb_var = mb_var, mc_var = mc_var, correlation = cor_val)
    })
  })
}

compute_r2 <- function(df, comparisons) {
  map_dfr(names(comparisons), function(mb_var) {
    map_dfr(comparisons[[mb_var]], function(mc_var) {
      fit <- tryCatch(lm(df[[mb_var]] ~ df[[mc_var]]), error = function(e) NULL)
      tibble(
        mb_var = mb_var,
        mc_var = mc_var,
        r2 = if (!is.null(fit)) summary(fit)$r.squared else NA_real_
      )
    })
  })
}

# ---- 4. Compute by Genus ----
cor_genus <- psmelt_all %>%
  group_by(Genus) %>%
  group_modify(~ compute_correlations(.x, comparisons = comparisons_log)) %>%
  ungroup()

r2_genus <- psmelt_all %>%
  group_by(Genus) %>%
  group_modify(~ compute_r2(.x, comparisons = comparisons_log)) %>%
  ungroup()

# ---- 5. Prepare plotting data ----
mc_mapping <- c(
  mc_counts_log = "Abundance", mc_relative_log = "Abundance",
  mc_carbon_log = "Carbon", carbon_relative_log = "Carbon",
  mc_volume_log = "Biovolume", biovol_relative_log = "Biovolume"
)

mb_labels <- c(
  mb_reads_log    = "Reads",
  mb_spike_log    = "Spike-norm",
  mb_dna_log      = "DNA-norm",
  mb_relative_log = "Relative (%)"
)

# Filter and map labels
cor_genus_plot <- cor_genus %>%
  filter(
    (mb_var %in% c("mb_reads_log", "mb_spike_log", "mb_dna_log") &
       mc_var %in% c("mc_counts_log", "mc_carbon_log", "mc_volume_log")) |
      (mb_var == "mb_relative_log" &
         mc_var %in% c("mc_relative_log", "carbon_relative_log", "biovol_relative_log"))
  ) %>%
  mutate(
    mc_group = mc_mapping[mc_var],
    mb_label = mb_labels[mb_var],
    mc_group = factor(mc_group, levels = c("Abundance", "Carbon", "Biovolume"))
  )

r2_genus_plot <- r2_genus %>%
  filter(
    (mb_var %in% c("mb_reads_log", "mb_spike_log", "mb_dna_log") &
       mc_var %in% c("mc_counts_log", "mc_carbon_log", "mc_volume_log")) |
      (mb_var == "mb_relative_log" &
         mc_var %in% c("mc_relative_log", "carbon_relative_log", "biovol_relative_log"))
  ) %>%
  mutate(
    mc_group = mc_mapping[mc_var],
    mb_label = mb_labels[mb_var],
    mc_group = factor(mc_group, levels = c("Abundance", "Carbon", "Biovolume"))
  )

# Ensure facet order and Genus ordering
panel_order <- c("Reads", "Spike-norm", "DNA-norm", "Relative (%)")

cor_genus_plot <- cor_genus_plot %>%
  mutate(
    mb_label = factor(mb_label, levels = panel_order),
    Genus = factor(Genus, levels = rev(sort(unique(Genus))))
  )

r2_genus_plot <- r2_genus_plot %>%
  mutate(
    mb_label = factor(mb_label, levels = panel_order),
    Genus = factor(Genus, levels = rev(sort(unique(Genus))))
  )

mean <- r2_genus_plot %>%
  group_by(mb_var, mc_var) %>%
  summarise(mean = mean(r2, na.rm = TRUE),
            .groups = "drop")

# Prepare the mean row
mean_plot <- mean %>%
  mutate(
    Genus = "Mean",          # new "species" row
    mc_group = case_when(    # map mc_var to mc_group
      grepl("counts", mc_var) ~ "Abundance",
      grepl("carbon", mc_var) ~ "Carbon",
      grepl("volume", mc_var) ~ "Biovolume",
      grepl("biovol", mc_var) ~ "Biovolume",
      grepl("relative", mc_var) ~ "Abundance",
      TRUE ~ "Other"
    ),
    mb_label = case_when(    # map mb_var to your labels
      mb_var == "mb_reads_log" ~ "Reads",
      mb_var == "mb_spike_log" ~ "Spike-norm",
      mb_var == "mb_dna_log" ~ "DNA-norm",
      mb_var == "mb_relative_log" ~ "Relative (%)",
      TRUE ~ "Other"
    ),
    r2 = mean                 # rename column for consistency
  ) %>%
  select(Genus, mb_var, mc_var, r2, mc_group, mb_label)

plot_data <- bind_rows(r2_genus_plot, mean_plot)

# Preserve original genus order
original_genera <- levels(r2_genus_plot$Genus)

# Factor levels: Mean last
plot_data$Genus <- factor(plot_data$Genus, levels = c(rev(original_genera), "Mean"))

# Reverse the levels for plotting so Mean is at the bottom
plot_data$Genus <- factor(plot_data$Genus, levels = rev(levels(plot_data$Genus)))

# ---- 6. Plotting ----
p1 <- ggplot(cor_genus_plot, aes(x = mc_group, y = Genus, fill = correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(correlation, 2)), size = 2.5, fontface = "bold") +
  facet_wrap(~mb_label, ncol = 4) +
  scale_fill_viridis(direction = -1, limits = c(0, 1), name = "Spearman correlation",
                     guide = guide_colorbar(direction = "horizontal",
                                            title.position = "top",
                                            title.hjust = 0.5)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 8, face = "bold.italic"),
    axis.text.x = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 11),
    legend.title.align = 0.5,
    legend.position = "bottom"
  ) +
  labs(title = "Genus-level Spearman correlations", x = "", y = "")

plot_data$mb_label <- factor(plot_data$mb_label,
                             levels = c("Reads", "Spike-norm", "DNA-norm", "Relative (%)"))

# prepare plotmath labels for y axis (keeps Mean non-italic)
genus_lvls <- levels(plot_data$Genus)
# escape single quotes if any names contain them
safe_lvls <- gsub("'", "\\\\'", genus_lvls)

labels_vec <- ifelse(
  genus_lvls == "Mean",
  paste0("bold('", safe_lvls, "')"),
  paste0("bold(italic('", safe_lvls, "'))")
)

p2 <- ggplot(plot_data, aes(x = mc_group, y = Genus, fill = r2)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(r2, 2)), size = 2.5, fontface = "bold") +
  facet_wrap(~mb_label, ncol = 4) +
  scale_y_discrete(labels = parse(text = labels_vec)) +           # <- use these labels
  scale_fill_viridis(direction = -1, limits = c(0, 1), name = expression(R^2),
                     guide = guide_colorbar(direction = "vertical",
                                            title.position = "top",
                                            title.hjust = 0.5)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 8),   # remove face here; plotmath handles style
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 11),
    legend.title.align = 0.5
  ) +
  labs(title = "", x = "", y = "")

ggsave("plots/microscopy_vs_barcoding/Supplementary_Figure_S3.png",
       p2,
       width = 8.27,   # in inches
       units = "in",
       dpi = 300)

