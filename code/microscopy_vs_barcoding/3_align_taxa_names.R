library(tidyverse)
library(vegan)

## Define paths
asv_seq_file = 'data/18S/ASVs_18S.txt' ## ASV sequences
seqtab_file = 'data/18S/seqtab_18S.txt' ## count matrix for each ASV & sample
microscopy_file = "data/shark_phytoplankton/phytoplankton_2019_2020_2023-06-02_utf8.txt"
# taxa_file = 'data/18S/taxa_18S.txt' ## taxonomic annotation for each ASV
taxa_file = 'data/18S/PR2_18S_2019_2020/raw/taxa_18S.txt' ## taxonomic annotation for each ASV
metadata_file = "data/18S/Balticdata_seqsamples_metadata_utf8.txt" ## metadata file (sample names, stations, dates etc.)
# translate_taxa_file = "data/18S/matching_taxa_microsc_vs_barcoding_17_march_2022.csv"
translate_taxa_file_sonia = "data/18S/MC_18S_nanomicrophyto_class_list_20230505.txt"
taxa_correction_file = 'data/18S/taxa_18S_PR2_5_1_0_SelectedClasses_SpCorrection_20250506.txt'

## Read preproccesing output
asv = read.table(asv_seq_file)$x
seqtab = as.matrix(read.delim(seqtab_file))
taxa = read.delim(taxa_file, sep = ' ', row.names = NULL)

# Read class list
translate_taxa_new = read.table(translate_taxa_file_sonia, fill = TRUE, header = TRUE, sep = "\t", na.strings="") %>%
  filter(MB_status == "included")

# Read list of taxa corrections
corrected_taxa = read.table(taxa_correction_file, fill = TRUE, header = TRUE, sep = "\t", na.strings="") %>%
  dplyr::select(#-ASV_nb,
    -AddSpecies,
    -Domain,
    -Supergroup,
    -Division,
    -Subdivision,
    -Class,
    -Order,
    -Family,
    -Genus,
    -Species,
    # -Comment
  ) %>%
  dplyr::rename(scientific_name = Scientific.name) %>%
  mutate(ASV = as.character(ASV)) %>%
  # mutate(ASV = paste0("ASV", ASV)) %>%
  mutate_all(na_if, "NA")

# Translate taxa table
taxa = taxa %>%
  left_join(corrected_taxa, by = c("row.names" = "ASV")) %>%
  dplyr::select(-Class, -Order, -Family, -Genus, -Species) %>%
  dplyr::rename(Class = CorrClass) %>%
  dplyr::rename(Order = CorrOrder) %>%
  dplyr::rename(Family = CorrFamily) %>%
  dplyr::rename(Genus = CorrGenus) %>%
  dplyr::rename(Species = CorrSpecies) %>%
  column_to_rownames(var = "row.names") %>%
  as.matrix()

## Read metadata
metadata = read.table(metadata_file, sep = '\t', header = TRUE)
meta_ix = match(colnames(seqtab), metadata$NGI_sample_ID_18S)
metadata = metadata[meta_ix,]

# remove blanks from the dataset
ix = sapply(
  'Blank', grep,
  metadata$station_name, invert = TRUE
)

# ix <- grep("Blank", metadata$station_name, invert = TRUE)


## Remove volume test

# NOTE: new meta_ix
meta_ix = intersect(which(metadata$origin != 'volume test'), ix)

## Merge stations B3 and B7 (practically one location)
metadata$station_name[metadata$station_name == 'NB1 / B3'] = 'B7'

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

## Choose only relevant samples for further use
metadata$station_date = paste(metadata$station_ID, metadata$sampling_date)

microscopy_data <- read.table(microscopy_file,
                              header = TRUE,
                              sep = "\t",
                              na.strings = "",
                              comment.char = "", # needed to avoid problems with "# counted"
                              fileEncoding = "utf8") # may need to be specified)

microscopy_data$station_date = paste(microscopy_data$station_id, microscopy_data$sample_date)


# Select only samples where we have microscopy data
meta_ix = intersect(which(metadata$station_date %in% microscopy_data$station_date), ix)


metadata = metadata[meta_ix,]

# Pick only selected stations for further diversity analyses
seqtab = seqtab[,meta_ix]
sample_id = metadata$sample_ID
colnames(seqtab) = sample_id

# Find zero reads and remove taxa
zero_reads = which(!rowSums(seqtab) == 0)

taxa = taxa[zero_reads,]
seqtab = seqtab[zero_reads,]
asv = asv[zero_reads]

### Get spike-in counts ###
## Identify spike ASVs

spike = 'CTTCGTTATCGTCACGGAGAGAACTGCCTTTAGCGATCTGTCTAGAGACCCGCCAAATATAAGCCTTGGGGTTATTAACAAAGACGCCAAACACTAGTGAATATGACAGTGAAGGGGTGTGAGATAGCTTTACTGGGTGAGAAAACACTCGTTAAAAAGAATTAGACCGGATAATCCCCGAGGGGCCGTAGGCATGGACTTGTCGTTGCCACCGAGCATAGCGGTTTCGAAATAGCCGAGATGGGCACTGGCGAATTAACCCACTGGTTTATATGGATCCGATGGGTTCACTTAATAAGCTCGTACCAGGGATGAATAAAGCGTTACGAGAATTATAAACATGGAGTTCCTATTGATTTGAGGTTAATACCGAACGGGAACATTTGTCGATCATGCTTCACATAGAGT'

spike_ix = agrep(spike, asv)

# ix = which(is.na(taxa[,2]))

taxa[spike_ix,]

spiketab = seqtab[spike_ix, ]

## Create dataframe with spike counts and rel data
## (including total reads and spike proportion)
rowSums(spiketab)
spike_counts = colSums(spiketab)
df_spike_counts = data.frame(sample_id = metadata$sample_ID_DNA[meta_ix],
                             station= metadata$station_name[meta_ix],
                             sampling_date = metadata$sampling_date[meta_ix],
                             # concentration = metadata$[meta_ix],
                             spike_cs = spike_counts)

tot_reads = colSums(seqtab)

df_spike_counts$tot_reads = tot_reads
df_spike_counts$spike_prop = df_spike_counts$spike_cs/tot_reads

### Chosing ASVs to inlude in further analyses ###
selected_classes = translate_taxa_new$Class_MB

## Exlude spike-ins from taxa annotation
ix_taxa =  setdiff(1:nrow(taxa), spike_ix)

#Excluding Metazoa
ix_taxa = intersect(ix_taxa, which(taxa[,4] != 'Metazoa'))
ix_taxa = intersect(ix_taxa, which(taxa[,1] != 'Eukaryota:nucl'))
ix_taxa = union(ix_taxa, which(is.na(taxa[,4])))

# Include only selected classes
ix_taxa = intersect(ix_taxa, which(taxa[,6] %in% selected_classes))

taxa = taxa[ix_taxa,]

taxa_all = as.data.frame(cbind(row.names(taxa), taxa)) %>%
  dplyr::rename("ASV" = V1) %>%
  mutate(Class = gsub("Prymnesiophyceae", "Coccolithophyceae", Class)) %>%
  mutate_all(funs(gsub(pattern = "_", replacement = " ", .)))

# Save table
write.table(taxa_all, "data/18S/taxa_18S_microscopy_new.txt", sep = "\t", row.names = FALSE)

