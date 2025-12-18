library(tidyverse)
library(phyloseq)
library(cowplot)
library(MicrobiotaProcess)
library(microbiome)
library(vegan)
library(pairwiseAdonis)

# File paths
phyloseq_file = 'data/18S/phyloseqs.RData'

# Load phyloseq objects
load(phyloseq_file)

# Transform to relative abundance (does not make a difference)
phyloseq_mb = transform_sample_counts(phyloseq_mb, function(x)100* x / sum(x))
phyloseq_abutab = transform_sample_counts(phyloseq_abutab, function(x)100* x / sum(x)) # Abundance
phyloseq_carbtab = transform_sample_counts(phyloseq_carbtab, function(x)100* x / sum(x)) # Carbon

# Select Abundance (phyloseq_abutab) or Carbon (phyloseq_carbtab)
phyloseq_um_selected <- phyloseq_carbtab

# Distmethod
pcoa_seqtab <- get_pcoa(obj=phyloseq_mb, distmethod="bray", method="hellinger")
pcoa_umtab <- get_pcoa(obj=phyloseq_um_selected, distmethod="bray", method="hellinger")

# Rename sea basins
pcoa_umtab@sampleda = pcoa_umtab@sampleda %>%
  mutate(sea_basin = dplyr::recode(sea_basin, 'Bothnian Bay' = 'I: Bothnian Bay', 'Bothnian Sea' = 'II: Bothnian Sea', 'Baltic Proper' = 'III: Baltic Proper', 'Kattegat' = 'IV: Kattegat', 'Skagerrak' = 'V: Skagerrak'))

pcoa_seqtab@sampleda = pcoa_seqtab@sampleda %>%
  mutate(sea_basin = dplyr::recode(sea_basin, 'Bothnian Bay' = 'I: Bothnian Bay', 'Bothnian Sea' = 'II: Bothnian Sea', 'Baltic Proper' = 'III: Baltic Proper', 'Kattegat' = 'IV: Kattegat', 'Skagerrak' = 'V: Skagerrak'))

# Visualize
pcoa_seqtab.p <- ggordpoint(obj=pcoa_seqtab, mapping = aes(alpha = 0.6), biplot=TRUE, speciesannot=TRUE,
                        factorNames=c("sea_basin", "season"), ellipse=TRUE, topn = 10) +
  labs(fill = "Sea basin", starshape = "Season") +
  ggtitle("Metabarcoding") +
  scale_alpha_continuous(guide = 'none')

# Extract ASV numbers
taxa = pcoa_seqtab.p[["layers"]][[4]][["data"]][["tax"]]

selected_taxa = as.data.frame(phyloseq::tax_table(phyloseq_mb)) %>%
  rownames_to_column() %>%
  filter(rowname %in% taxa) %>%
  arrange(rowname, levels = taxa)

selected_taxa = selected_taxa[match(taxa, selected_taxa$rowname), ]

# Replace ASV numbers with taxa names
pcoa_seqtab.p[["layers"]][[4]][["data"]][["tax"]] = paste(taxa, selected_taxa$scientific_name, sep = ": ")

pcoa_umtab.p <- ggordpoint(obj=pcoa_umtab, mapping = aes(alpha = 0.6), biplot=TRUE, speciesannot=TRUE,
                        factorNames=c("sea_basin", "season"), ellipse=TRUE, topn = 10) +
  labs(fill = "Sea basin", starshape = "Season") +
  ggtitle("Utermöhl") +
  scale_alpha_continuous(guide = 'none')

# Extract plot ranges
x_min = min(c(ggplot_build(pcoa_seqtab.p)$layout$panel_scales_x[[1]]$range$range,
               ggplot_build(pcoa_umtab.p)$layout$panel_scales_x[[1]]$range$range))

x_max = max(c(ggplot_build(pcoa_seqtab.p)$layout$panel_scales_x[[1]]$range$range,
              ggplot_build(pcoa_umtab.p)$layout$panel_scales_x[[1]]$range$range))

y_min = min(c(ggplot_build(pcoa_seqtab.p)$layout$panel_scales_y[[1]]$range$range,
                ggplot_build(pcoa_umtab.p)$layout$panel_scales_y[[1]]$range$range))

y_max = max(c(ggplot_build(pcoa_seqtab.p)$layout$panel_scales_y[[1]]$range$range,
              ggplot_build(pcoa_umtab.p)$layout$panel_scales_y[[1]]$range$range))

x_range = c(x_min, x_max)
y_range = c(y_min, y_max)

# Set limits for both plots
pcoa_umtab.p <- pcoa_umtab.p +
  xlim(x_range) +
  ylim(y_range)

pcoa_seqtab.p <- pcoa_seqtab.p +
  xlim(x_range) +
  ylim(y_range)

# Extract the legend from one of the plots
legend <- ggfun::get_legend(
  pcoa_seqtab.p + guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          legend.margin = margin(t = -0.5, unit='cm'),
          legend.title = element_blank())
)

# Plots in grid
plot = plot_grid(pcoa_seqtab.p + theme(legend.position = "none"),
                 pcoa_umtab.p + theme(legend.position = "none"),
                 legend,
                 align = "v",
                 axis = "l",
                 labels = c("a)", "b)"),
                 label_x = 0.08,
                 label_y = 0.9,
                 ncol = 1,
                 nrow = 3,
                 rel_heights =  c(1, 1, .2))

# Save plot
ggsave("fig7_pcoa_carbon.png",
       width = 7,
       height = 7,
       plot = plot,
       device = "png",
       path = "plots/microscopy_vs_barcoding/",
       bg = "white")


plot_ordination(
  physeq = phyloseq_mb,
  ordination = pcoa_seqtab,
  type="scree")

# PERMANOVA

# MB
metadata <- as(sample_data(phyloseq_mb), "data.frame")

permanova.mb = adonis2(distance(phyloseq_mb, method="bray") ~ sea_basin,
       data = metadata)

beta <- betadisper(distance(phyloseq_mb, method="bray"), metadata$sea_basin)
permutest(beta)
boxplot_mb <- boxplot(beta, xlab=NA, main = "Metabarcoding")

pairwise.adonis(distance(phyloseq_mb, method="bray"), metadata$sea_basin)

# SIMPER
# simper_mb <- simper(t(otu_table(phyloseq_mb)), sample_data(phyloseq_mb)$sea_basin)

write.table(summary(simper_mb)$Kattegat_Skagerrak, "data/tables/simper_mb.txt", sep = "\t", row.names = TRUE)

# ANOSIM
mb_basin <- get_variable(phyloseq_mb, "sea_basin")

anosim_mb = anosim(distance(phyloseq_mb, "bray"), mb_basin)

# UM
metadata_um <- as(sample_data(phyloseq_um_selected), "data.frame")

permanova.um = adonis2(distance(phyloseq_um_selected, method="bray") ~ sea_basin,
                       data = metadata_um)

beta_um <- betadisper(distance(phyloseq_um_selected, method="bray"), metadata_um$sea_basin)
permutest(beta_um)
boxplot_um <- boxplot(beta_um, xlab=NA, main = "Utermöhl")

pairwise.adonis(distance(phyloseq_um_selected, method="bray"), metadata_um$sea_basin)

# SIMPER
# simper_um <- simper(t(otu_table(phyloseq_um_selected)), sample_data(phyloseq_um_selected)$sea_basin)

write.table(summary(simper_um)$Kattegat_Skagerrak, "data/tables/simper_um.txt", sep = "\t", row.names = TRUE)

# ANOSIM
um_basin <- get_variable(phyloseq_um_selected, "sea_basin")

anosim_um = anosim(distance(phyloseq_um_selected, "bray"), um_basin)
