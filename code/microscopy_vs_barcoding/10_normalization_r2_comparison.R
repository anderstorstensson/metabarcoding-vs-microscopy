library(tidyverse)
library(phyloseq)
library(cowplot)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(purrr)
library(patchwork)
library(viridis)

# File paths
phyloseq_file = 'data/18S/phyloseqs.RData'

# Load phyloseq objects
load(phyloseq_file)

# Subset samples to identical samples
phyloseq_abutab = prune_samples(sample_data(phyloseq_abutab)$NGI_sample_ID_18S %in%
                                  sample_data(phyloseq_spike)$NGI_sample_ID_18S, phyloseq_abutab)

phyloseq_carbtab = prune_samples(sample_data(phyloseq_carbtab)$NGI_sample_ID_18S %in%
                                  sample_data(phyloseq_spike)$NGI_sample_ID_18S, phyloseq_carbtab)

phyloseq_voltab = prune_samples(sample_data(phyloseq_voltab)$NGI_sample_ID_18S %in%
                                   sample_data(phyloseq_spike)$NGI_sample_ID_18S, phyloseq_voltab)

phyloseq_mb = prune_samples(sample_data(phyloseq_mb)$NGI_sample_ID_18S %in%
                                  sample_data(phyloseq_spike)$NGI_sample_ID_18S, phyloseq_mb)

phyloseq_dna = prune_samples(sample_data(phyloseq_dna)$NGI_sample_ID_18S %in%
                              sample_data(phyloseq_spike)$NGI_sample_ID_18S, phyloseq_dna)

# Aggregate TO class level
phyloseq_mb = tax_glom(phyloseq_mb, taxrank="Class", NArm = FALSE)
phyloseq_spike = tax_glom(phyloseq_spike, taxrank="Class", NArm = FALSE)
phyloseq_dna = tax_glom(phyloseq_dna, taxrank="Class", NArm = FALSE)
phyloseq_abutab = tax_glom(phyloseq_abutab, taxrank="Class_MB", NArm = FALSE)
phyloseq_carbtab = tax_glom(phyloseq_carbtab, taxrank="Class_MB", NArm = FALSE)
phyloseq_voltab = tax_glom(phyloseq_voltab, taxrank="Class_MB", NArm = FALSE)

# Transform to relative abundance
phyloseq_mb_rel = transform_sample_counts(phyloseq_mb, function(x)100* x / sum(x))
phyloseq_abutab_rel = transform_sample_counts(phyloseq_abutab, function(x)100* x / sum(x))
phyloseq_carbtab_rel = transform_sample_counts(phyloseq_carbtab, function(x)100* x / sum(x))
phyloseq_voltab_rel = transform_sample_counts(phyloseq_voltab, function(x)100* x / sum(x))

# Melt tables
psmelt_seqtab = psmelt(phyloseq_mb) %>%
  filter(!Abundance == 0) %>%
  mutate(method = "Metabarcoding") %>%
  mutate(abundance_type = "Reads") %>%
  select(NGI_sample_ID_18S,
         sea_basin,
         Class,
         Abundance) %>%
  rename("mb_reads" = Abundance)

psmelt_seqtab_rel = psmelt(phyloseq_mb_rel) %>%
  filter(!Abundance == 0) %>%
  mutate(method = "Metabarcoding") %>%
  mutate(abundance_type = "Relative abundance") %>%
  select(NGI_sample_ID_18S,
         sea_basin,
         Class,
         Abundance) %>%
  rename("mb_relative" = Abundance)

psmelt_spiketab = psmelt(phyloseq_spike) %>%
  filter(!Abundance == 0) %>%
  mutate(method = "Metabarcoding") %>%
  mutate(abundance_type = "Spike-normalized") %>%
  select(NGI_sample_ID_18S,
         sea_basin,
         Class,
         Abundance) %>%
  rename("mb_spike" = Abundance)

psmelt_dnatab = psmelt(phyloseq_dna) %>%
  filter(!Abundance == 0) %>%
  mutate(method = "Metabarcoding") %>%
  mutate(abundance_type = "DNA-normalized") %>%
  select(NGI_sample_ID_18S,
         sea_basin,
         Class,
         Abundance) %>%
  rename("mb_dna" = Abundance)

psmelt_abutab = psmelt(phyloseq_abutab) %>%
  filter(!Abundance == 0) %>%
  mutate(method = "Microscopy") %>%
  mutate(abundance_type = "Cell counts") %>%
  dplyr::rename(Class = Class_MB) %>%
  select(NGI_sample_ID_18S,
         sea_basin,
         Class,
         Abundance) %>%
  rename("mc_counts" = Abundance)

psmelt_carbtab = psmelt(phyloseq_carbtab) %>%
  filter(!Abundance == 0) %>%
  mutate(method = "Microscopy") %>%
  mutate(abundance_type = "Carbon concentration") %>%
  dplyr::rename(Class = Class_MB) %>%
  select(NGI_sample_ID_18S,
         sea_basin,
         Class,
         Abundance) %>%
  rename("mc_carbon" = Abundance)

psmelt_voltab = psmelt(phyloseq_voltab) %>%
  filter(!Abundance == 0) %>%
  mutate(method = "Microscopy") %>%
  mutate(abundance_type = "Biovolume concentration") %>%
  dplyr::rename(Class = Class_MB) %>%
  select(NGI_sample_ID_18S,
         sea_basin,
         Class,
         Abundance) %>%
  rename("mc_volume" = Abundance)

psmelt_abutab_rel = psmelt(phyloseq_abutab_rel) %>%
  filter(!Abundance == 0) %>%
  mutate(method = "Microscopy") %>%
  mutate(abundance_type = "Relative abundance") %>%
  dplyr::rename(Class = Class_MB) %>%
  select(NGI_sample_ID_18S,
         sea_basin,
         Class,
         Abundance) %>%
  rename("mc_relative" = Abundance)

psmelt_carbtab_rel = psmelt(phyloseq_carbtab_rel) %>%
  filter(!Abundance == 0) %>%
  mutate(method = "Microscopy") %>%
  mutate(abundance_type = "Carbon relative abundance") %>%
  dplyr::rename(Class = Class_MB) %>%
  select(NGI_sample_ID_18S,
         sea_basin,
         Class,
         Abundance) %>%
  rename("carbon_relative" = Abundance)

psmelt_voltab_rel = psmelt(phyloseq_voltab_rel) %>%
  filter(!Abundance == 0) %>%
  mutate(method = "Microscopy") %>%
  mutate(abundance_type = "Biovolume relative abundance") %>%
  dplyr::rename(Class = Class_MB) %>%
  select(NGI_sample_ID_18S,
         sea_basin,
         Class,
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

psmelt_combined = psmelt_all %>%
  mutate(sea_basin = "All")

psmelt_all = rbind(psmelt_all, psmelt_combined)

# Define comparisons as a list
comparisons <- list(
  mb_reads  = c("mc_counts", "mc_carbon", "mc_volume"),
  mb_spike  = c("mc_counts", "mc_carbon", "mc_volume"),
  mb_dna    = c("mc_counts", "mc_carbon", "mc_volume"),
  mb_relative = c("mc_relative", "carbon_relative", "biovol_relative")
)

# Function to compute correlations
compute_correlations <- function(df) {
  map_dfr(names(comparisons), function(mb_var) {
    map_dfr(comparisons[[mb_var]], function(mc_var) {
      cor_val <- cor(df[[mb_var]], df[[mc_var]], use = "complete.obs", method = "spearman")
      tibble(mb_var = mb_var, mc_var = mc_var, correlation = cor_val)
    })
  })
}

# --- General correlations (all data)
cor_general <- compute_correlations(psmelt_all)

# --- By sea basin
cor_basin <- psmelt_all %>%
  group_by(sea_basin) %>%
  group_modify(~ compute_correlations(.x)) %>%
  ungroup()

# --- By class
cor_class <- psmelt_all %>%
  group_by(Class) %>%
  group_modify(~ compute_correlations(.x)) %>%
  ungroup()



# Define the full desired order for both axes
mb_levels <- c("mb_reads", "mb_spike", "mb_dna", "mb_relative")
mc_levels <- c("mc_counts", "mc_carbon", "mc_volume",
               "mc_relative", "carbon_relative", "biovol_relative")

# Mapping for nicer labels and ordering
mb_labels <- c(
  mb_reads    = "Reads",
  mb_spike    = "Spike-norm",
  mb_dna      = "DNA-norm",
  mb_relative = "Relative (%)"
)

# Prepare the plotting dataframe
cor_plot <- cor_general %>%
  filter(
    (mb_var %in% c("mb_reads", "mb_spike", "mb_dna") & mc_var %in% c("mc_counts", "mc_carbon", "mc_volume")) |
      (mb_var == "mb_relative" & mc_var %in% c("mc_relative", "carbon_relative", "biovol_relative"))
  ) %>%
  mutate(
    # Group absolute + relative metrics into columns
    mc_group = case_when(
      mc_var %in% c("mc_counts", "mc_relative") ~ "Counts",
      mc_var %in% c("mc_carbon", "carbon_relative") ~ "Carbon",
      mc_var %in% c("mc_volume", "biovol_relative") ~ "Biovolume"
    ),
    # Rename and order mb_var
    mb_label = factor(mb_labels[mb_var], levels = mb_labels[mb_levels]),
    mc_group = factor(mc_group, levels = c("Counts", "Carbon", "Biovolume"))
  )

cor_plot <- cor_plot %>%
  mutate(
    mb_label = factor(mb_label, levels = rev(levels(mb_label)))  # Reverse for ggplot
  )

# Plot
ggplot(cor_plot, aes(x = mc_group, y = mb_label, fill = correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(correlation, 2)), size = 3) +
  scale_fill_viridis(
    direction = -1,
    limits = c(0.4, 0.75), name = "Spearman\ncorrelation"
  ) +
  theme_minimal(base_size = 13) +
  labs(
    title = "General correlations between metabarcoding and microscopy parameters",
    x = "Microscopy group",
    y = "Metabarcoding variable"
  )

# --- Function to compute R² values
compute_r2 <- function(df) {
  map_dfr(names(comparisons), function(mb_var) {
    map_dfr(comparisons[[mb_var]], function(mc_var) {
      # Fit linear model: mb_var ~ mc_var
      fit <- lm(df[[mb_var]] ~ df[[mc_var]])
      tibble(
        mb_var = mb_var,
        mc_var = mc_var,
        r2 = summary(fit)$r.squared
      )
    })
  })
}

# --- Compute general R²
r2_general <- compute_r2(psmelt_all)

# --- Prepare the plotting dataframe for R2 (reuse same mapping & ordering)
r2_plot <- r2_general %>%
  filter(
    (mb_var %in% c("mb_reads", "mb_spike", "mb_dna") & mc_var %in% c("mc_counts", "mc_carbon", "mc_volume")) |
      (mb_var == "mb_relative" & mc_var %in% c("mc_relative", "carbon_relative", "biovol_relative"))
  ) %>%
  mutate(
    mc_group = case_when(
      mc_var %in% c("mc_counts", "mc_relative") ~ "Counts",
      mc_var %in% c("mc_carbon", "carbon_relative") ~ "Carbon",
      mc_var %in% c("mc_volume", "biovol_relative") ~ "Biovolume"
    ),
    mb_label = factor(mb_labels[mb_var], levels = mb_labels[mb_levels]),
    mc_group = factor(mc_group, levels = c("Counts", "Carbon", "Biovolume"))
  ) %>%
  mutate(
    mb_label = factor(mb_label, levels = rev(levels(mb_label)))  # Reverse for ggplot
  )



# --- Spearman plot
p1 <- ggplot(cor_plot, aes(x = mc_group, y = mb_label, fill = correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(correlation, 2)), size = 4, fontface = "bold") +
  scale_fill_viridis(
    direction = -1,
    limits = c(0.4, 0.75),
    name = "Spearman correlation",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", angle = 0, hjust = 0.5),
    axis.text.y = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title.align = 0.5,
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(title = "Spearman correlations", x = NULL, y = NULL)

# --- R² plot
p2 <- ggplot(r2_plot, aes(x = mc_group, y = mb_label, fill = r2)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(r2, 2)), size = 4, fontface = "bold") +
  scale_fill_viridis(
    direction = -1,
    limits = c(0, 0.75),
    name = expression(R^2),
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", angle = 0, hjust = 0.5),
    axis.text.y = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title.align = 0.5,
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(title = "Linear regression R²", x = NULL, y = NULL)

# --- Combine side by side
combined <- p1 + p2 + plot_layout(ncol = 2)

ggsave("plots/microscopy_vs_barcoding/Fig8.png", combined)



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

# ---- 2. Define log-transformed comparisons ----
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

# ---- 4. Compute correlations and R² ----
cor_general <- compute_correlations(psmelt_all, comparisons_log)
r2_general  <- compute_r2(psmelt_all, comparisons_log)

# ---- 5. Define labels and ordering ----
mb_labels <- c(
  mb_reads_log    = "Reads",
  mb_spike_log    = "Spike-norm",
  mb_dna_log      = "DNA-norm",
  mb_relative_log = "Relative (%)"
)

mc_mapping <- c(
  mc_counts_log = "Abundance", mc_relative_log = "Abundance",
  mc_carbon_log = "Carbon", carbon_relative_log = "Carbon",
  mc_volume_log = "Biovolume", biovol_relative_log = "Biovolume"
)

# ---- 6. Prepare plotting dataframes ----
prepare_plot_df <- function(df) {
  df %>%
    filter(
      (mb_var %in% c("mb_reads_log", "mb_spike_log", "mb_dna_log") &
         mc_var %in% c("mc_counts_log", "mc_carbon_log", "mc_volume_log")) |
        (mb_var == "mb_relative_log" &
           mc_var %in% c("mc_relative_log", "carbon_relative_log", "biovol_relative_log"))
    ) %>%
    mutate(
      mc_group = mc_mapping[mc_var],
      mb_label = factor(mb_labels[mb_var], levels = mb_labels),
      mc_group = factor(mc_group, levels = c("Abundance", "Carbon", "Biovolume"))
    ) %>%
    mutate(
      mb_label = factor(mb_label, levels = rev(levels(mb_label)))  # Reverse for ggplot
    )
}

cor_plot <- prepare_plot_df(cor_general)
r2_plot  <- prepare_plot_df(r2_general)

# ---- 7. Plotting ----
# Spearman correlation plot
p1 <- ggplot(cor_plot, aes(x = mc_group, y = mb_label, fill = correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(correlation, 2)), size = 4, fontface = "bold") +
  scale_fill_viridis(
    direction = -1,
    limits = c(0.4, 0.75),
    name = "Spearman correlation",
    guide = guide_colorbar(direction = "horizontal", title.position = "top", title.hjust = 0.5)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", angle = 0, hjust = 0.5),
    axis.text.y = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title.align = 0.5,
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(title = "Spearman correlations", x = NULL, y = NULL)

# R² plot
p2 <- ggplot(r2_plot, aes(x = mc_group, y = mb_label, fill = r2)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(r2, 2)), size = 4, fontface = "bold") +
  scale_fill_viridis(
    direction = -1,
    limits = c(0, 0.75),
    name = expression(R^2),
    guide = guide_colorbar(direction = "horizontal", title.position = "top", title.hjust = 0.5)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", angle = 0, hjust = 0.5),
    axis.text.y = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title.align = 0.5,
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(title = "Linear regression R²", x = NULL, y = NULL)

# Combine plots side by side
combined <- p1 + p2 + plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a")  # adds "a" and "b" labels automatically

# Print the combined plot
combined

ggsave("plots/microscopy_vs_barcoding/Fig8_log.png", combined, width = 8.27, height = 4)
