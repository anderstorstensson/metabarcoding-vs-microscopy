library(tidyverse)
library(cowplot)
library(phyloseq)

# File path
phyloseq_file = 'data/18S/phyloseqs.RData' 

# Load phyloseq objects
load(phyloseq_file)

# Get taxa lists for MB and UM
taxa_all = as.data.frame(tax_table(phyloseq_mb)) %>%
  rownames_to_column(var = "ASV")

taxa_um = as.data.frame(tax_table(phyloseq_abutab))
taxa_um[taxa_um == ''] <- NA

# Prepare MB tables
barcoding_table_match = taxa_all %>%
  dplyr::group_by(Class) %>%
  dplyr::summarise(Taxaincommon.species = sum((na.omit(unique(Species)) %in% unique(taxa_um$taxon_species))==TRUE),
                   Taxaincommon.genus = sum((na.omit(unique(Genus)) %in% unique(taxa_um$taxon_genus))==TRUE),
                   Taxaincommon.order = sum((na.omit(unique(Order)) %in% unique(taxa_um$taxon_order))==TRUE))


barcoding_table = taxa_all %>%
  dplyr::group_by(Class) %>%
  dplyr::summarise(ASVs = length(ASV),
                   ASVs_species = sum(!is.na(Species), na.rm=TRUE),
                   ASVs_order = sum(!is.na(Order), na.rm=TRUE),
                   ASVs_genus = sum(!is.na(Genus), na.rm=TRUE),
                   NumberofSpecies = sum(!is.na(unique(Species))),
                   NumberofOrder = sum(!is.na(unique(Order))),
                   NumberofGenera = sum(!is.na(unique(Genus)))) %>%
  left_join(barcoding_table_match)

# Prepare UM tables
microscopy_table = taxa_um %>%
  dplyr::group_by(Class_MB) %>%
  dplyr::summarise(NumberofSpecies = sum(!is.na(unique(taxon_species))),
                   NumberofOrder = sum(!is.na(unique(taxon_order))),
                   NumberofGenera = sum(!is.na(unique(taxon_genus))))

# Join tables
full_table = microscopy_table %>%
  dplyr::rename("Class" = "Class_MB") %>%
  right_join(barcoding_table, by = "Class")

# Cleanup of tables
for(i in 1:ncol(full_table)) {
  names(full_table)[i] = gsub("\\.x", ".microscopy", names(full_table)[i])
}

for(i in 1:ncol(full_table)) {
  names(full_table)[i] = gsub("\\.y", ".barcoding", names(full_table)[i])
}

full_table = full_table %>%
  arrange(Class) %>%
  replace(is.na(.), 0)

# Save table
write.table(full_table, file = "data/18S/comparison_table.txt", row.names = FALSE, sep = "\t", na = "")

# Prepare tables for plotting
all_in_one_class = full_table %>%
  # filter(rank == "Class") %>%
  gather("Method", "Value",-Class) %>%
  dplyr::rename("Level" = Class) %>%
  mutate("Parameter" = Method) %>%
  mutate(Method = gsub("NumberofOrder.microscopy", "Microscopy", Method)) %>%
  mutate(Method = gsub("NumberofGenera.microscopy", "Microscopy", Method)) %>%
  mutate(Method = gsub("NumberofSpecies.microscopy", "Microscopy", Method)) %>%
  mutate(Method = gsub("NumberofOrder.barcoding", "Metabarcoding", Method)) %>%
  mutate(Method = gsub("NumberofGenera.barcoding", "Metabarcoding", Method)) %>%
  mutate(Method = gsub("NumberofSpecies.barcoding", "Metabarcoding", Method)) %>%
  mutate(Method = gsub("Taxaincommon.species", "Taxa in common", Method)) %>%
  mutate(Method = gsub("Taxaincommon.genus", "Taxa in common", Method)) %>%
  mutate(Method = gsub("Taxaincommon.order", "Taxa in common", Method)) %>%
  mutate(Parameter = gsub("NumberofOrder.microscopy", "c) Order", Parameter)) %>%
  mutate(Parameter = gsub("NumberofGenera.microscopy", "b) Genus", Parameter)) %>%
  mutate(Parameter = gsub("NumberofSpecies.microscopy", "a) Species", Parameter)) %>%
  mutate(Parameter = gsub("NumberofOrder.barcoding", "c) Order", Parameter)) %>%
  mutate(Parameter = gsub("NumberofGenera.barcoding", "b) Genus", Parameter)) %>%
  mutate(Parameter = gsub("NumberofSpecies.barcoding", "a) Species", Parameter)) %>%
  mutate(Parameter = gsub("Taxaincommon.species", "a) Species", Parameter)) %>%
  mutate(Parameter = gsub("Taxaincommon.genus", "b) Genus", Parameter)) %>%
  mutate(Parameter = gsub("Taxaincommon.order", "c) Order", Parameter)) %>%
  filter(Method == "Microscopy" | Method == "Metabarcoding" | Method == "Taxa in common") %>%
  mutate(rank = "Class")

write.table(all_in_one_class, file = "data/18S/comparison_table_all.txt", row.names = FALSE, sep = "\t", na = "")

# Wrangle dataframes for plot
taxa_keep = all_in_one_class %>%
  select(Level, Value) %>%
  mutate(Value = as.numeric(Value)) %>%
  dplyr::group_by(Level) %>%
  dplyr::summarise(ntaxa = sum(Value, na.rm = TRUE)) %>%
  arrange(ntaxa) 

all_in_one_class = all_in_one_class %>%
  filter(Level %in% taxa_keep$Level) %>%
  mutate(Level = as.factor(Level)) %>%
  mutate(Level = gsub("ophyceae", "", Level)) %>%
  mutate(Method = gsub("Microscopy", "UM only", Method)) %>%
  mutate(Method = gsub("Metabarcoding", "MB only", Method)) %>%
  mutate(Method = gsub("Taxa in common", "MB & UM", Method)) %>%
  mutate(Method = factor(Method, levels = c("MB only", "UM only", "MB & UM")))
  
all_in_one_class = all_in_one_class %>%
  mutate(Level = gsub("Coscinodisc + Medi", "Coscinodisc\n+ Medi", Level, fixed = TRUE))

plot = ggplot(all_in_one_class, aes(x = fct_reorder(Level, desc(as.numeric(Value))), y = as.integer(Value), fill = Method)) +
  geom_bar(position="dodge", 
           stat="identity", 
           color = "black",
           linewidth = 0.25) +
  theme_bw() +
  theme(
        legend.position="top",
        legend.title = element_blank(),
        # panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        text = element_text(size = 10),
        legend.key.height = unit(.1, 'cm'),
        strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("Number of identified taxa") +
  xlab("") +
  scale_fill_manual(values=c('white', '#999999', 'black')) +
  facet_grid(cols = vars(factor(rank, levels = c("Phylum", "Class"))), 
             rows = vars(factor(Parameter, levels = c("a) Species", "b) Genus", "c) Order", "Taxa in common"))), 
             scales = "free")
  # scale_x_discrete(guide = guide_axis(n.dodge = 2))

ggsave("fig3_number_of_taxa.png", 
       width = 5,
       height = 5,
       plot = plot, 
       device = "png", 
       path = "plots/microscopy_vs_barcoding/")
