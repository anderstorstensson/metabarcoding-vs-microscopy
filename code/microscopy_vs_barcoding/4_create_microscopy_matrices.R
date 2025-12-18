library(tidyverse)
library(readxl)

## Define paths

seqtab_file = 'data/18S/seqtab_18S.txt' ## count matrix for each ASV & sample
microscopy_file = "data/shark_phytoplankton/phytoplankton_2019_2020_2023-06-02_utf8.txt"
metadata_file = "data/18S/Balticdata_seqsamples_metadata.txt" ## metadata file (sample names, stations, dates etc.)
translate_taxa_file_sonia = "data/18S/MC_18S_nanomicrophyto_class_list_20230505.txt"
translate_taxa_file_anders = "data/18S/translations.txt"
taxa_correction_file_microscopy_sonia = 'data/18S/Taxa list 20250204.xlsx'

translations_anders <- read_tsv(translate_taxa_file_anders)
microscopy_corrections_sonia <- read_excel(taxa_correction_file_microscopy_sonia) %>%
  mutate(scientific_name_corrected = scientific_name) %>%
  select(-taxon_class,
         -taxon_order,
         -taxon_family,
         -taxon_genus,
         -trophic_type_code,
         -bvol_scientific_name,
         -bvol_aphia_id,
         -scientific_name,
         -dyntaxa_id,
         -aphia_id,
         -worms_class,
         -worms_order,
         -worms_family,
         -worms_genus) %>%
  distinct()

microscopy_data <- read.table(microscopy_file,
                              header = TRUE,
                              sep = "\t",
                              na.strings = "",
                              comment.char = "", # needed to avoid problems with "# counted"
                              fileEncoding = "utf8") # may need to be specified)

microscopy_data <- microscopy_data %>%
  left_join(microscopy_corrections_sonia) %>%
  select(-scientific_name,
         -taxon_phylum,
         -taxon_class,
         -taxon_order,
         -taxon_family,
         -taxon_genus) %>%
  rename(scientific_name = scientific_name_corrected,
         taxon_phylum = ab_phylum,
         taxon_class = ab_class,
         taxon_order = ab_order,
         taxon_family = ab_family,
         taxon_genus = ab_genus)

microscopy_data <- microscopy_data %>%
  mutate(taxon_genus = gsub("Plagioselmis", "Teleaulax", taxon_genus)) %>%
  mutate(scientific_name = gsub("Plagioselmis prolonga", "Teleaulax amphioxeia", scientific_name)) %>%
  mutate(scientific_name = gsub("Plagioselmis", "Teleaulax", scientific_name)) %>%
  mutate(taxon_species = gsub("Plagioselmis prolonga", "Teleaulax amphioxeia", taxon_species)) %>%
  mutate(taxon_species = gsub("Emiliania huxleyi", "Gephyrocapsa huxleyi", taxon_species)) %>%
  mutate(scientific_name = gsub("Chaetoceros ceratosporus var. ceratosporus", "Chaetoceros ceratosporus", scientific_name)) %>%
  mutate(scientific_name = gsub("Chaetoceros subtilis var. subtilis", "Chaetoceros subtilis", scientific_name)) %>%
  mutate(scientific_name = gsub("Chaetoceros throndsenii var. throndsenii", "Chaetoceros throndsenii", scientific_name)) %>%
  mutate(scientific_name = gsub("Tripos horridus", "Tripos longipes", scientific_name)) %>%
  mutate(taxon_species = gsub("Tripos horridus", "Tripos longipes", taxon_species)) %>%
  mutate(scientific_name = gsub("Volvocales", "Chlamydomonadales", scientific_name)) %>%
  mutate(taxon_class = if_else(
    taxon_class %in% c("Coscinodiscophyceae", "Mediophyceae"),
    "Coscinodiscophyceae + Mediophyceae",
    taxon_class
  )) %>%
  left_join(translations_anders, by = "scientific_name") %>%
  mutate(scientific_name = if_else(!is.na(translate_to), translate_to, scientific_name)) %>%
  mutate(taxon_species = if_else(!is.na(species) & species, translate_to, taxon_species)) %>%
  select(-species, -translate_to) %>%  # Remove the temporary translation column
  mutate(scientific_name = if_else(reported_scientific_name == "Chaetoceros gracilis", "Chaetoceros gracilis", scientific_name)) %>%
  mutate(taxon_species = if_else(reported_scientific_name == "Chaetoceros gracilis", "Chaetoceros gracilis", taxon_species))

translate_taxa_new <- read.table(translate_taxa_file_sonia, fill = TRUE, header = TRUE, sep = "\t", na.strings="") %>%
  filter(MB_status == "included")

translate <- translate_taxa_new %>%
  dplyr::select(Class_MB, Group_MC)

# Summarise number of unidentified unicells and flagellates
uni_flag = microscopy_data %>%
  mutate(station_date = paste(station_id, sample_date)) %>%
  filter(scientific_name == "Unicell" | scientific_name == "Flagellates") %>%
  filter(parameter == "Abundance") %>%
  dplyr::select(station_date,
                value,
                scientific_name) %>%
  dplyr::group_by(station_date, scientific_name) %>%
  dplyr::summarise(value = sum(value)) %>%
  pivot_wider(names_from = scientific_name, values_from = value)

microscopy_data <- microscopy_data %>%
  left_join(translate, by = c("taxon_class" = "Group_MC")) %>%
  left_join(translate, by = c("taxon_phylum" = "Group_MC")) %>%
  mutate("Class_MB" = coalesce(Class_MB.x, Class_MB.y)) %>%
  dplyr::select(-Class_MB.x, -Class_MB.y) %>%
  filter(!is.na(Class_MB)) %>%
  mutate(Class_MB = gsub("Prymnesiophyceae", "Coccolithophyceae", Class_MB)) %>%
  mutate(taxon_class = Class_MB)

## Read preproccesing output

seqtab <- as.matrix(read.delim(seqtab_file))

metadata <- read.table(metadata_file, sep = '\t', header = TRUE, fileEncoding = "macroman")
meta_ix <-  match(colnames(seqtab), metadata$NGI_sample_ID_18S)

metadata <-  metadata[meta_ix,]

# remove blanks from the dataset

ix <- sapply(
  'Blank', grep,
  metadata$station_name, invert = TRUE
)

## Remove volume test

# NOTE: new meta_ix
meta_ix <- intersect(which(metadata$origin != 'volume test'), ix)

## Merge stations B3 and B7 (practically one location)

metadata$station_name[metadata$station_name == 'NB1 / B3'] = 'B7'

# Remove oversequenced outlier
meta_ix <- meta_ix[meta_ix != 2]
meta_ix <- meta_ix[meta_ix != 134]

# Choose just one of replicates from the same time and station
dates_stations <- c()

for(i in meta_ix){
  ds <- paste(metadata$sampling_date[i], metadata$station_name[i])
  if(! ds %in% dates_stations){
    dates_stations <- c(dates_stations, ds)
  }else{
    meta_ix <- meta_ix[! meta_ix %in% i]
  }
}

## Choose only relevant metadate for further use
metadata <- metadata[meta_ix,]

metadata$station_date <- paste(metadata$station_ID, metadata$sampling_date)

microscopy_data$station_date <- paste(microscopy_data$station_id, microscopy_data$sample_date)

# Select only relevant samples
microscopy_data <- microscopy_data %>%
  filter(station_date %in% metadata$station_date)

# Check samples without match in the two datasets
no_microscopy_match <- metadata %>%
  filter(!metadata$station_date %in% microscopy_data$station_date)

# Aggregate size class info
abundance_data <- microscopy_data %>%
  filter(parameter == "Abundance") %>%
  dplyr::select(-size_class,
                -magnification,
                -coefficient,
                -reported_value) %>%
  dplyr::group_by(across(c(-value))) %>%
  dplyr::summarise(value = sum(value))

# Aggregate size class info
biovolume_data <- microscopy_data %>%
  filter(parameter == "Biovolume concentration") %>%
  dplyr::select(-size_class,
                -magnification,
                -coefficient,
                -reported_value) %>%
  dplyr::group_by(across(c(-value))) %>%
  dplyr::summarise(value = sum(value))

# Aggregate size class info
carbon_data <- microscopy_data %>%
  filter(parameter == "Carbon concentration") %>%
  dplyr::select(-size_class,
                -magnification,
                -coefficient,
                -reported_value) %>%
  dplyr::group_by(across(c(-value))) %>%
  dplyr::summarise(value = sum(value))

# Aggregate size class info
count_data <- microscopy_data %>%
  filter(parameter == "# counted") %>%
  dplyr::select(-size_class,
                -magnification,
                -coefficient,
                -bvol_size_class,
                -reported_value) %>%
  dplyr::group_by(across(c(-value))) %>%
  dplyr::summarise(value = sum(value))

# Transform abundance data into matrix
species_matrix <- abundance_data %>%
  ungroup() %>%
  dplyr::select(scientific_name, station_date, value) %>%
  pivot_wider(
    names_from = station_date,
    values_from = value,
    values_fn = sum,
    values_fill = 0
  )

species_matrix <- as.data.frame(species_matrix)

row.names(species_matrix) <- species_matrix$scientific_name

species_matrix <- species_matrix %>%
  select(-scientific_name)

names <- data.frame(names(species_matrix)) %>%
  dplyr::rename("station_date" = names.species_matrix.) %>%
  left_join(metadata, by = "station_date")

names(species_matrix) <- names$NGI_sample_ID_18S

# Transform biovolume data into matrix
biovolume_matrix <- biovolume_data %>%
  ungroup() %>%
  dplyr::select(scientific_name, station_date, value) %>%
  pivot_wider(
    names_from = station_date,
    values_from = value,
    values_fn = sum,
    values_fill = 0
  )

biovolume_matrix <- as.data.frame(biovolume_matrix)

row.names(biovolume_matrix) <- biovolume_matrix$scientific_name

biovolume_matrix <- biovolume_matrix %>%
  select(-scientific_name)

names <- data.frame(names(biovolume_matrix)) %>%
  dplyr::rename("station_date" = names.biovolume_matrix.) %>%
  left_join(metadata, by = "station_date")

names(biovolume_matrix) <- names$NGI_sample_ID_18S

# Transform carbon data into matrix
carbon_matrix <- carbon_data %>%
  ungroup() %>%
  dplyr::select(scientific_name, station_date, value) %>%
  pivot_wider(
    names_from = station_date,
    values_from = value,
    values_fn = sum,
    values_fill = 0
  )

carbon_matrix <- as.data.frame(carbon_matrix)

row.names(carbon_matrix) <- carbon_matrix$scientific_name

carbon_matrix <- carbon_matrix %>%
  select(-scientific_name)

names <- data.frame(names(carbon_matrix)) %>%
  dplyr::rename("station_date" = names.carbon_matrix.) %>%
  left_join(metadata, by = "station_date")

names(carbon_matrix) <- names$NGI_sample_ID_18S

# Transform count data into matrix
count_matrix <- count_data %>%
  ungroup() %>%
  dplyr::select(scientific_name, station_date, value) %>%
  pivot_wider(
    names_from = station_date,
    values_from = value,
    values_fn = sum,
    values_fill = 0
  )

count_matrix <- as.data.frame(count_matrix)

row.names(count_matrix) <- count_matrix$scientific_name

count_matrix <- count_matrix %>%
  select(-scientific_name)

names <- data.frame(names(count_matrix)) %>%
  dplyr::rename("station_date" = names.count_matrix.) %>%
  left_join(metadata, by = "station_date")

names(count_matrix) <- names$NGI_sample_ID_18S

# Create a taxa table
tax_table <- abundance_data %>%
  ungroup() %>%
  dplyr::select(Class_MB,
                taxon_kingdom,
         taxon_phylum,
         taxon_class,
         taxon_order,
         taxon_family,
         taxon_genus,
         taxon_species,
         scientific_name) %>%
  distinct()

tax_table <- as.data.frame(tax_table)

row.names(tax_table) <- tax_table$scientific_name

tax_table <- tax_table %>%
  select(-scientific_name)

# Add unicells and flagellates to metadata
names <- names %>%
  left_join(uni_flag)

# Save matrices
write.table(count_matrix, file = "data/shark_phytoplankton/microscopy_count.txt", row.names = TRUE, sep = "\t", na = "")
write.table(species_matrix, file = "data/shark_phytoplankton/microscopy_abundance.txt", row.names = TRUE, sep = "\t", na = "")
write.table(biovolume_matrix, file = "data/shark_phytoplankton/microscopy_biovolume.txt", row.names = TRUE, sep = "\t", na = "")
write.table(carbon_matrix, file = "data/shark_phytoplankton/microscopy_carbon.txt", row.names = TRUE, sep = "\t", na = "")
write.table(tax_table, file = "data/shark_phytoplankton/tax_table.txt", row.names = TRUE, sep = "\t", na = "")
write.table(names, file = "data/shark_phytoplankton/microscopy_metadata.txt", row.names = TRUE, sep = "\t", na = "")
