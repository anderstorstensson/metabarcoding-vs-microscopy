library(tidyverse)
library(vegan)
library(phyloseq)

# File paths
asv_seq_file = 'data/18S/ASVs_18S.txt' ## ASV sequences
seqtab_file = 'data/18S/seqtab_18S.txt' ## count matrix for each ASV & sample
taxa_file = 'data/18S/taxa_18S.txt' ## taxonomic annotation for each ASV
metadata_file = "data/18S/Balticdata_seqsamples_metadata.txt" ## metadata file (sample names, stations, dates etc.)
metadata_file_microscopy = "data/shark_phytoplankton/microscopy_metadata.txt" ## metadata file (sample names, stations, dates etc.)
meta_ix_file = "data/18S/meta_ix.csv" # indices of samples that we had sequencing data for
normalized_seqtab_file = "data/18S/seqtab_18S_spike_normalised.txt"
spike_metadata_file = "data/18S/spike_normalized_metadata.txt"
abundance_file =  'data/shark_phytoplankton/microscopy_abundance.txt'
biovolume_file =  'data/shark_phytoplankton/microscopy_biovolume.txt'
carbon_file =  'data/shark_phytoplankton/microscopy_carbon.txt'
taxonomy_file = 'data/shark_phytoplankton/tax_table.txt'
translate_taxa_file_sonia = "data/18S/MC_18S_nanomicrophyto_class_list_20230505.txt"
read_abundance_file = "data/18S/read_abundance.txt"
taxa_correction_file = 'data/18S/taxa_18S_PR2_5_01_AssignTaxonomy_SpCorrection_20230607.txt'

# Read translate
translate_taxa_new = read.table(translate_taxa_file_sonia, fill = TRUE, header = TRUE, sep = "\t", na.strings="") %>%
  filter(MB_status == "included")

## Read preproccesing output
asv = read.table(asv_seq_file)$x
abutab = as.matrix(read.delim(abundance_file))
voltab = as.matrix(read.delim(biovolume_file))
carbtab = as.matrix(read.delim(carbon_file))
taxa = as.matrix(read.delim(taxonomy_file, sep = '\t'))
seqtab = as.matrix(read.delim(seqtab_file))
spiketab = as.matrix(read.delim(normalized_seqtab_file))
taxa_mb = read.delim(taxa_file, sep = ' ', row.names = NULL)

# Read taxa translation file
corrected_taxa = read.table(taxa_correction_file, fill = TRUE, header = TRUE, sep = "\t", na.strings="") %>%
  dplyr::select(-ASV_nb, 
                -Domain, 
                -Supergroup, 
                -Division,
                -Subdivision, 
                -Class,
                -Order,
                -Family, 
                -Genus,
                -Species,
                -Comment) %>%
  dplyr::rename(scientific_name = Scientific.name) %>%
  mutate_all(na_if, "NA")

# Correct taxonomy
taxa_mb = taxa_mb %>%
  left_join(corrected_taxa, by = c("row.names" = "ASV")) %>%
  dplyr::select(-Class, -Order, -Family, -Genus, -Species) %>%
  dplyr::rename(Class = CorrClass) %>%
  dplyr::rename(Order = CorrOrder) %>%
  dplyr::rename(Family = CorrFamily) %>%
  dplyr::rename(Genus = CorrGenus) %>%
  dplyr::rename(Species = CorrSpecies) %>%
  column_to_rownames(var = "row.names") %>%
  as.matrix()

# Read microscopy metadata
metadata_microscopy = read.table(metadata_file_microscopy, sep = '\t', header = TRUE)

# Read read abundance file
read_abundance = read.table(read_abundance_file, sep = '\t', header = TRUE)

# Read metabarcoding metadata
metadata = read.table(metadata_file, sep = '\t', header = TRUE) %>%
  left_join(read_abundance)
metadata_spike = read.table(spike_metadata_file, sep = '\t', header = TRUE)

# Match metadata and sequence abundance table 
meta_ix = match(colnames(seqtab), metadata$NGI_sample_ID_18S)
meta_ix_spike = match(colnames(spiketab), metadata_spike$NGI_sample_ID_18S)

# Select only metadata that are present in the sequence abundance table
metadata = metadata[meta_ix,]
metadata_spike = metadata_spike[meta_ix_spike,]

# remove blanks from the datasets
ix = sapply(
  'Blank', grep,
  metadata$station_name, invert = TRUE
)

ix_spike = sapply(
  'Blank', grep,
  metadata_spike$station_name, invert = TRUE
)

## Remove volume test

# NOTE: new meta_ix
meta_ix = intersect(which(metadata$origin != 'volume test'), ix)
meta_ix_spike = intersect(which(metadata_spike$origin != 'volume test'), ix)

## Merge stations B3 and B7 (practically one location)
metadata$station_name[metadata$station_name == 'NB1 / B3'] = 'B7'
metadata_spike$station_name[metadata_spike$station_name == 'NB1 / B3'] = 'B7'

# Remove oversequenced outlier
meta_ix = meta_ix[meta_ix != 2]
meta_ix = meta_ix[meta_ix != 134]

# Choose just one of replicates from the same time and station
dates_stations = c()

for(i in meta_ix){
  ds = paste(metadata$sampling_date[i], metadata$station_name[i])
  if(! ds %in% dates_stations){
    dates_stations = c(dates_stations, ds)
  }else{
    meta_ix = meta_ix[! meta_ix %in% i]
  }
}

dates_stations_spike = c()

for(i in meta_ix){
  ds = paste(metadata_spike$sampling_date[i], metadata_spike$station_name[i])
  if(! ds %in% dates_stations_spike){
    dates_stations_spike = c(dates_stations_spike, ds)
  }else{
    meta_ix_spike = meta_ix_spike[! meta_ix_spike %in% i]
  }
}

# Match metadata and sequence abundance table 
meta_ix = intersect(which(metadata$NGI_sample_ID_18S %in% metadata_microscopy$NGI_sample_ID_18S), ix)
meta_ix_spike = intersect(which(metadata_spike$NGI_sample_ID_18S %in% metadata_microscopy$NGI_sample_ID_18S), ix_spike)

# Select only metadata that are present in the sequence abundance table
metadata = metadata[meta_ix,]
metadata_spike = metadata_spike[meta_ix_spike,]

# Pick only selected samples for further analyses
seqtab = seqtab[,meta_ix]
spiketab = spiketab[,meta_ix_spike]

## Identify spike ASVs
spike = 'CTTCGTTATCGTCACGGAGAGAACTGCCTTTAGCGATCTGTCTAGAGACCCGCCAAATATAAGCCTTGGGGTTATTAACAAAGACGCCAAACACTAGTGAATATGACAGTGAAGGGGTGTGAGATAGCTTTACTGGGTGAGAAAACACTCGTTAAAAAGAATTAGACCGGATAATCCCCGAGGGGCCGTAGGCATGGACTTGTCGTTGCCACCGAGCATAGCGGTTTCGAAATAGCCGAGATGGGCACTGGCGAATTAACCCACTGGTTTATATGGATCCGATGGGTTCACTTAATAAGCTCGTACCAGGGATGAATAAAGCGTTACGAGAATTATAAACATGGAGTTCCTATTGATTTGAGGTTAATACCGAACGGGAACATTTGTCGATCATGCTTCACATAGAGT'

# Perform approximate string matching between the spike sequence and another set of sequences called ASVs (amplicon sequence variants).
# The result is stored in the variable 'spike_ix'.
spike_ix = agrep(spike, asv)

## Exlude spike-ins and Metazoa from taxa annotation

ix_taxa =  setdiff(1:nrow(taxa_mb), spike_ix)
ix_taxa = intersect(ix_taxa, which(taxa_mb[,4] != 'Metazoa'))
ix_taxa = intersect(ix_taxa, which(taxa_mb[,1] != 'Eukaryota:nucl'))

# Subset tables
taxa_mb = taxa_mb[ix_taxa,]
seqtab = seqtab[ix_taxa,]
spiketab = spiketab[ix_taxa,]

# Rarefying seqtab
m = min(colSums(seqtab))
m
r_seqtab = rrarefy(x = t(seqtab), sample = m)
r_seqtab = t(r_seqtab)

# Chosing ASVs to inlude in further analyses
selected_classes = translate_taxa_new$Class_MB
ix_taxa =  1:nrow(taxa_mb)

# Include only selected classes
ix_taxa = intersect(ix_taxa, which(taxa_mb[,5] %in% selected_classes))

# Subset tables
taxa_mb = taxa_mb[ix_taxa,]
seqtab = seqtab[ix_taxa,]
r_seqtab = r_seqtab[ix_taxa,]
spiketab = spiketab[ix_taxa,]

# Find zero reads
zero_reads = which(!rowSums(seqtab) == 0)
zero_reads_r_seqtab = which(!rowSums(r_seqtab) == 0)
zero_read_spiketab = which(!rowSums(spiketab) == 0)

# Remove zero read taxa in tax_table
taxa_spike = taxa_mb[zero_read_spiketab,]
taxa_mb_r = taxa_mb[zero_reads_r_seqtab,]
taxa_mb = taxa_mb[zero_reads,]

# Remove zero read taxa in seqtabs
seqtab = seqtab[zero_reads,]
r_seqtab = r_seqtab[zero_reads_r_seqtab,]
spiketab = spiketab[zero_read_spiketab,]

# Create metadata table compatible with PS
metadata_ps = metadata
rownames(metadata_ps) = metadata_ps$NGI_sample_ID_18S

# Translate station name to sea basin
metadata_ps = metadata_ps %>%
  mutate(sea_basin = station_name_alt)

# Define a lookup table for replacements
replacement_table = data.frame(
  pattern = c("RA2", "RA1", "A13", "B3", "B7", "NB1 / B3", "C3", "GA1", "SR3/C24", "B1",
              "BY31", "BY15", "Ref M1V1", "BY2", "BY5", "BCS III-10", "BY38", "Anholt E",
              "N14 Falkenberg", "Sl\x8agg\x9a", "A17", "\x8117"),
  replace_with = c("Bothnian Bay", "Bothnian Bay", "Bothnian Bay", "Bothnian Sea",
                   "Bothnian Sea", "Bothnian Sea", "Bothnian Sea", "Bothnian Sea",
                   "Bothnian Sea", "Baltic Proper", "Baltic Proper", "Baltic Proper",
                   "Baltic Proper", "Baltic Proper", "Baltic Proper", "Baltic Proper",
                   "Baltic Proper", "Kattegat", "Kattegat", "Skagerrak",
                   "Skagerrak", "Skagerrak")
)

# Loop through the lookup table and perform replacements
for (i in 1:nrow(replacement_table)) {
  metadata_ps$sea_basin = gsub(replacement_table$pattern[i], replacement_table$replace_with[i],
                                metadata_ps$sea_basin, ignore.case = FALSE, perl = FALSE,
                                fixed = FALSE, useBytes = FALSE)
}

# Transform into factor
metadata_ps$sea_basin = factor(metadata_ps$sea_basin, levels = c("Bothnian Bay", 
                                                                 "Bothnian Sea", 
                                                                 "Baltic Proper", 
                                                                 "Kattegat", 
                                                                 "Skagerrak"))
# Add day of year to metadata
metadata_ps$day_of_year = yday(metadata_ps$sampling_date)
metadata_ps$day_of_year = factor(metadata_ps$day_of_year)
metadata_ps$day_of_year = fct_reorder(metadata_ps$day_of_year, as.integer(metadata_ps$day_of_year))

# Add season to metadata
metadata_ps = metadata_ps %>%
  mutate(season = case_when(
    sampling_month %in% c("1", "2", "12") ~ "Dec-Feb",
    sampling_month %in% c("3", "4", "5") ~ "Mar-May",
    sampling_month %in% c("6", "7", "8") ~ "Jun-Aug",
    sampling_month %in% c("9", "10", "11") ~ "Sep-Nov"
  )) %>%
  mutate(season = factor(season, levels = c("Dec-Feb", "Mar-May", "Jun-Aug", "Sep-Nov")))

# Add method to metadata
metadata_ps$method = "Metabarcoding"

# Do the same for spike data

# Create metadata table compatible with PS
metadata_ps_spike = metadata_spike
rownames(metadata_ps_spike) = metadata_ps_spike$NGI_sample_ID_18S

# Translate station name to sea basin
metadata_ps_spike = metadata_ps_spike %>%
  mutate(sea_basin = station_name_alt)

# Define a lookup table for replacements
replacement_table = data.frame(
  pattern = c("RA2", "RA1", "A13", "B3", "B7", "NB1 / B3", "C3", "GA1", "SR3/C24", "B1",
              "BY31", "BY15", "Ref M1V1", "BY2", "BY5", "BCS III-10", "BY38", "Anholt E",
              "N14 Falkenberg", "Sl\xe4gg\xf6", "A17", "\xc517"),
  replace_with = c("Bothnian Bay", "Bothnian Bay", "Bothnian Bay", "Bothnian Sea",
                   "Bothnian Sea", "Bothnian Sea", "Bothnian Sea", "Bothnian Sea",
                   "Bothnian Sea", "Baltic Proper", "Baltic Proper", "Baltic Proper",
                   "Baltic Proper", "Baltic Proper", "Baltic Proper", "Baltic Proper",
                   "Baltic Proper", "Kattegat", "Kattegat", "Skagerrak",
                   "Skagerrak", "Skagerrak")
)

# Loop through the lookup table and perform replacements
for (i in 1:nrow(replacement_table)) {
  metadata_ps_spike$sea_basin = gsub(replacement_table$pattern[i], replacement_table$replace_with[i],
                                      metadata_ps_spike$sea_basin, ignore.case = FALSE, perl = FALSE,
                                      fixed = FALSE, useBytes = FALSE)
}

# Transform into factor
metadata_ps_spike$sea_basin = factor(metadata_ps_spike$sea_basin, levels = c("Bothnian Bay", 
                                                                 "Bothnian Sea", 
                                                                 "Baltic Proper", 
                                                                 "Kattegat", 
                                                                 "Skagerrak"))
# Add day of year to metadata
metadata_ps_spike$day_of_year = yday(metadata_ps_spike$sampling_date)
metadata_ps_spike$day_of_year = factor(metadata_ps_spike$day_of_year)
metadata_ps_spike$day_of_year = fct_reorder(metadata_ps_spike$day_of_year, as.integer(metadata_ps_spike$day_of_year))

# Add season to metadata
metadata_ps_spike = metadata_ps_spike %>%
  mutate(season = case_when(
    sampling_month %in% c("1", "2", "12") ~ "Dec-Feb",
    sampling_month %in% c("3", "4", "5") ~ "Mar-May",
    sampling_month %in% c("6", "7", "8") ~ "Jun-Aug",
    sampling_month %in% c("9", "10", "11") ~ "Sep-Nov"
  )) %>%
  mutate(season = factor(season, levels = c("Dec-Feb", "Mar-May", "Jun-Aug", "Sep-Nov")))

# Add method to metadata
metadata_ps_spike$method = "Metabarcoding"

# Add scientific name to taxa table and clean up names for plots
taxa_clean = as.data.frame(taxa_mb) %>%
  mutate_all(funs(gsub(pattern = "_", replacement = " ", .))) %>%
  as.matrix()

# Add scientific name to rarified taxa table and clean up names for plots
taxa_r_clean = as.data.frame(taxa_mb_r) %>%
  mutate_all(funs(gsub(pattern = "_", replacement = " ", .))) %>%
  as.matrix()

# Add scientific name to spike taxa table and clean up names for plots
taxa_spike_clean = as.data.frame(taxa_spike) %>%
  mutate_all(funs(gsub(pattern = "_", replacement = " ", .))) %>%
  as.matrix()

# Create DNA-normalized table to total read abundance (including all classes)

# Calculate relative abundance
norm_seqtab = seqtab
for (i in 1:ncol(seqtab)) {
  norm_seqtab[,i] = seqtab[,i]/metadata_ps$reads.nospikenometazoa[i]
}

# Multiply relative abundance with DNA concentration
seqtab_dna_norm = norm_seqtab %*% diag(metadata_ps$sample_DNA_concentration)

# Rename column names
colnames(seqtab_dna_norm) = metadata_ps$NGI_sample_ID_18S

# Create metabarcoding PS objects
phyloseq_mb = phyloseq(otu_table(seqtab, taxa_are_rows=TRUE), 
                     sample_data(metadata_ps), 
                     tax_table(as.matrix(taxa_clean))
)

phyloseq_spike = phyloseq(otu_table(spiketab, taxa_are_rows=TRUE), 
                        sample_data(metadata_ps_spike), 
                        tax_table(as.matrix(taxa_spike_clean))
)

phyloseq_dna = phyloseq(otu_table(seqtab_dna_norm, taxa_are_rows=TRUE), 
                           sample_data(metadata_ps), 
                           tax_table(as.matrix(taxa_clean))
)

phyloseq_r_mb = phyloseq(otu_table(r_seqtab, taxa_are_rows=TRUE), 
                         sample_data(metadata_ps), 
                         tax_table(as.matrix(taxa_r_clean))
)

# Create metadata table compatible with PS object from microscopy
metadata_ps = metadata_microscopy
rownames(metadata_ps) = metadata_ps$NGI_sample_ID_18S

# Add sea basin based on station name
metadata_ps = metadata_ps %>%
  mutate(sea_basin = station_name_alt)

# Define a lookup table for replacements
replacement_table = data.frame(
  pattern = c("RA2", "RA1", "A13", "B3", "B7", "NB1 / B3", "C3", "GA1", "SR3/C24", "B1",
              "BY31", "BY15", "Ref M1V1", "BY2", "BY5", "BCS III-10", "BY38", "Anholt E",
              "N14 Falkenberg", "Sl\x8agg\x9a", "SlŠggš", "A17", "\x8117"),
  replace_with = c("Bothnian Bay", "Bothnian Bay", "Bothnian Bay", "Bothnian Sea",
                   "Bothnian Sea", "Bothnian Sea", "Bothnian Sea", "Bothnian Sea",
                   "Bothnian Sea", "Baltic Proper", "Baltic Proper", "Baltic Proper",
                   "Baltic Proper", "Baltic Proper", "Baltic Proper", "Baltic Proper",
                   "Baltic Proper", "Kattegat", "Kattegat", "Skagerrak", "Skagerrak",
                   "Skagerrak", "Skagerrak")
)

# Loop through the lookup table and perform replacements
for (i in 1:nrow(replacement_table)) {
  metadata_ps$sea_basin = gsub(replacement_table$pattern[i], replacement_table$replace_with[i],
                                metadata_ps$sea_basin, ignore.case = FALSE, perl = FALSE,
                                fixed = FALSE, useBytes = FALSE)
}

# Transform to factor
metadata_ps$sea_basin = factor(metadata_ps$sea_basin, levels = c("Bothnian Bay", "Bothnian Sea", "Baltic Proper", "Kattegat", "Skagerrak"))

# Add day of year to metadata
metadata_ps$day_of_year = yday(metadata_ps$sampling_date)
metadata_ps$day_of_year = factor(metadata_ps$day_of_year)
metadata_ps$day_of_year = fct_reorder(metadata_ps$day_of_year, as.integer(metadata_ps$day_of_year))

# Add season to metadata
metadata_ps = metadata_ps %>%
  mutate(season = case_when(
    sampling_month %in% c("1", "2", "12") ~ "Dec-Feb",
    sampling_month %in% c("3", "4", "5") ~ "Mar-May",
    sampling_month %in% c("6", "7", "8") ~ "Jun-Aug",
    sampling_month %in% c("9", "10", "11") ~ "Sep-Nov"
  )) %>%
  mutate(season = factor(season, levels = c("Dec-Feb", "Mar-May", "Jun-Aug", "Sep-Nov")))

# Add method to metdata
metadata_ps$method = "Microscopy"


# Create PS object
phyloseq_abutab = phyloseq(otu_table(abutab, taxa_are_rows=TRUE), 
                     sample_data(metadata_ps), 
                     tax_table(as.matrix(taxa))
)

phyloseq_carbtab = phyloseq(otu_table(carbtab, taxa_are_rows=TRUE), 
                            sample_data(metadata_ps), 
                            tax_table(as.matrix(taxa))
)

phyloseq_voltab = phyloseq(otu_table(voltab, taxa_are_rows=TRUE), 
                             sample_data(metadata_ps), 
                             tax_table(as.matrix(taxa))
)

# Save all objects
save(phyloseq_mb, 
     phyloseq_spike, 
     phyloseq_r_mb, 
     phyloseq_dna, 
     phyloseq_abutab, 
     phyloseq_carbtab, 
     phyloseq_voltab, 
     file = "data/18S/phyloseqs.RData")
