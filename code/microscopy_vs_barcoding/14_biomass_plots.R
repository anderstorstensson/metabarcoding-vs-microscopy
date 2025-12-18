library(tidyverse)
library(phyloseq)

# File paths 
phyloseq_file = 'data/18S/phyloseqs.RData' 
mc_no_translate_file = 'data/shark_phytoplankton/phytoplankton_2019_2020_2023-06-17_utf8_no_taxa_translation.txt'
taxa_file = 'data/18S/taxa_18S.txt' ## taxonomic annotation for each ASV


### Metabarcoding table ###

# Load phyloseq objects
load(phyloseq_file)

cell_counts = psmelt(phyloseq_abutab) %>%
  group_by(NGI_sample_ID_18S) %>%
  summarise(Abundance = sum(Abundance)) %>%
  left_join(sample_data(phyloseq_abutab)) %>%
  mutate("parameter" = "Cell count")

biovolume = psmelt(phyloseq_voltab) %>%
  group_by(NGI_sample_ID_18S) %>%
  summarise(Abundance = sum(Abundance)) %>%
  left_join(sample_data(phyloseq_abutab)) %>%
  mutate("parameter" = "Biovolume")

carbon = psmelt(phyloseq_carbtab) %>%
  group_by(NGI_sample_ID_18S) %>%
  summarise(Abundance = sum(Abundance)) %>%
  left_join(sample_data(phyloseq_abutab)) %>%
  mutate("parameter" = "Carbon")

dna = as.tibble(sample_data(phyloseq_abutab)) %>%
  mutate("Abundance" = sample_DNA_concentration) %>%
  mutate("parameter" = "DNA")

all = rbind(cell_counts, biovolume, carbon, dna)

all %>%
  ggplot(aes(x = as.Date(sampling_date), y = Abundance, color = parameter)) +
  geom_point() +
  facet_wrap(vars(parameter), scales = "free_y") +
  scale_x_date(date_labels = "%b")
