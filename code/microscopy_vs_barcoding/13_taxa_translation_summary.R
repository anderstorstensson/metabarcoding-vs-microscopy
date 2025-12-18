library(tidyverse)
library(phyloseq)

# File paths 
phyloseq_file = 'data/18S/phyloseqs.RData' 
mc_no_translate_file = 'data/shark_phytoplankton/phytoplankton_2019_2020_2023-06-17_utf8_no_taxa_translation.txt'
taxa_file = 'data/18S/taxa_18S.txt' ## taxonomic annotation for each ASV


### Metabarcoding table ###

# Load phyloseq objects
load(phyloseq_file)

# Extract MB taxa table from phyloseq
taxa_mb_translated = as.data.frame(tax_table(phyloseq_mb)) %>%
  rownames_to_column(var = "ASV") %>%
  select(-scientific_name, -Domain, -Supergroup, -Division, -Subdivision)
taxa_mb_translated[taxa_mb_translated == ''] <- NA

# Read MB taxa file
taxa_mb_no_translate = read.delim(taxa_file, sep = ' ', row.names = NULL) %>%
  rename("ASV" = row.names) %>%
  mutate_all(funs(gsub(pattern = "_", replacement = " ", .))) %>%
  filter(ASV %in% taxa_mb_translated$ASV) %>%
  select(-Domain, -Supergroup, -Division, -Subdivision)

# Identify difference between the two taxa tables
difference = setdiff(taxa_mb_translated, taxa_mb_no_translate)

# Extract the taxa that differed
taxa_mb_no_translate = taxa_mb_no_translate %>%
  filter(ASV %in% difference$ASV) %>%
  rename("Class_from" = Class, "Order_from" = Order, "Family_from" = Family, "Genus_from" = Genus, "Species_from" = Species)

taxa_mb_translated = taxa_mb_translated %>%
  filter(ASV %in% difference$ASV) %>%
  rename("Class_to" = Class, "Order_to" = Order, "Family_to" = Family, "Genus_to" = Genus, "Species_to" = Species)

# Join tables
mb_translate_summary = taxa_mb_translated %>%
  left_join(taxa_mb_no_translate)

# Extract class translations
class = mb_translate_summary %>%
  select(Class_from, Class_to) %>%
  distinct() %>%
  filter(!Class_to == Class_from) %>%
  rename("To" = Class_to, "From" = Class_from) %>%
  mutate(Rank = "Class") %>%
  arrange(From)

# Extract order translations
order = mb_translate_summary %>%
  select(Order_from, Order_to) %>%
  distinct() %>%
  filter(!Order_to == Order_from) %>%
  rename("To" = Order_to, "From" = Order_from) %>%
  mutate(Rank = "Order") %>%
  arrange(From)

# Extract family translations
family = mb_translate_summary %>%
  select(Family_from, Family_to) %>%
  distinct() %>%
  filter(!Family_to == Family_from) %>%
  rename("To" = Family_to, "From" = Family_from) %>%
  mutate(Rank = "Family") %>%
  arrange(From)

# Extract genus translations
genus = mb_translate_summary %>%
  select(Genus_from, Genus_to) %>%
  distinct() %>%
  filter(!Genus_to == Genus_from) %>%
  rename("To" = Genus_to, "From" = Genus_from) %>%
  mutate(Rank = "Genus") %>%
  arrange(From)

# Extract species translations
species = mb_translate_summary %>%
  select(Species_from, Species_to) %>%
  distinct() %>%
  filter(!Species_to == Species_from) %>%
  rename("To" = Species_to, "From" = Species_from) %>%
  mutate(Rank = "Species") %>%
  arrange(From)

# Bind tables together
summary = rbind(class, order, family, genus, species)

# Save table
write.table(summary, "data/tables/mb_translation_summary.txt", sep = "\t", row.names = FALSE)

### Microscopy table ###

# Read non-translated microscopy file
microscopy_data_no_translate <- read.table(mc_no_translate_file,
                                           header = TRUE,
                                           sep = "\t",
                                           na.strings = "",
                                           comment.char = "", # needed to avoid problems with "# counted"
                                           fileEncoding = "utf8") # may need to be specified)



# Extract UM taxa table from phyloseq
taxa_translated = as.data.frame(tax_table(phyloseq_abutab)) %>%
  distinct()
taxa_translated[taxa_translated == ''] <- NA


# Summarise non-translated taxa list
taxa_non_translated = microscopy_data_no_translate %>%
  select(taxon_kingdom,
         taxon_phylum,
         taxon_class,
         taxon_order,
         taxon_family,
         taxon_genus,
         taxon_species) %>%
  distinct()

# Find changes
species_to = unique(taxa_translated$taxon_species[!taxa_translated$taxon_species %in% taxa_non_translated$taxon_species])
genus_to = unique(taxa_translated$taxon_genus[!taxa_translated$taxon_genus %in% taxa_non_translated$taxon_genus])
family_to = unique(taxa_translated$taxon_family[!taxa_translated$taxon_family %in% taxa_non_translated$taxon_family])
order_to = unique(taxa_translated$taxon_order[!taxa_translated$taxon_order %in% taxa_non_translated$taxon_order])
class_to = unique(taxa_translated$taxon_class[!taxa_translated$taxon_class %in% taxa_non_translated$taxon_class])

