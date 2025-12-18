library(tidyverse)
library(phyloseq)

## Define infiles
metadata_file = "data/18S/Balticdata_seqsamples_metadata.txt" ## metadata file (sample names, stations, dates etc.)
asv_seq_file = 'data/18S/ASVs_18S.txt' ## ASV sequences
seqtab_file = 'data/18S/seqtab_18S.txt' ## count matrix for each ASV & sample
taxa_file = 'data/18S/taxa_18S.txt' ## taxonomic annotation for each ASV
old_metadata_file = "data/18S/P20107_samples_metadata.csv"
filtration_volume_smhi_file = "data/18S/filtration_info_smhi.txt"

# Read filtration info
filtration_volumes = read.table(filtration_volume_smhi_file, fill = TRUE, header = TRUE, sep = "\t") %>%
  dplyr::rename("old_sample_id" = sample_ID_DNA,
                "filtration_volume" = FILTVOL) %>%
  dplyr::select(old_sample_id, filtration_volume) %>%
  distinct()

# Read and process old metadata file, renaming columns and joining with filtration volumes
metadata_old = read.csv(old_metadata_file, sep = ";") %>%
  dplyr::rename("NGI_sample_ID_18S" = NGI_sample_ID,
                "old_sample_id" = sample_ID) %>%
  dplyr::select(NGI_sample_ID_18S, old_sample_id) %>%
  mutate("old_sample_id" = gsub("_", " ", old_sample_id)) %>%
  left_join(filtration_volumes)

# Read current metadata file and join with old metadata
metadata = read.csv(metadata_file, sep = '\t') %>%
  left_join(metadata_old) 

# Loop through each row in metadata and calculate spike per liter based on origin
for(i in 1:nrow(metadata)) {
  if(metadata$origin[i] %in% c("SMHI", "volume test")) {
    metadata$spike_per_liter_18S[i] = 0.01/metadata$filtration_volume[i]
  }
  if(metadata$origin[i] == "UMU") {
    metadata$spike_per_liter_18S[i] = 0.007/metadata$filtration_volume[i]
  }
}

# Read ASV, sequence count matrix, and taxonomic annotation data
asv = read.table(asv_seq_file)$x
seqtab = as.matrix(read.delim(seqtab_file))
taxa = as.matrix(read.delim(taxa_file, sep = ' '))

# Filter metadata to include only samples present in the sequence count matrix
metadata = metadata %>%
  filter(NGI_sample_ID_18S %in% colnames(seqtab))

# Create an index of metadata rows
meta_ix = seq(1, nrow(metadata))

# Remove blanks from the dataset
ix = sapply(
  'Blank', grep,
  metadata$station_name, invert = TRUE
)

meta_ix = intersect(which(metadata$spike_per_liter_18S > 0), ix)
# meta_ix = intersect(meta_ix, which(metadata$spike_per_liter_18S > 0))

# Remove over-sequenced outliers
meta_ix = meta_ix[meta_ix != 2]
meta_ix = meta_ix[meta_ix != 134]

### Get spike-in counts ###

## Identify spike ASVs

# A DNA sequence representing the spike.
spike = 'CTTCGTTATCGTCACGGAGAGAACTGCCTTTAGCGATCTGTCTAGAGACCCGCCAAATATAAGCCTTGGGGTTATTAACAAAGACGCCAAACACTAGTGAATATGACAGTGAAGGGGTGTGAGATAGCTTTACTGGGTGAGAAAACACTCGTTAAAAAGAATTAGACCGGATAATCCCCGAGGGGCCGTAGGCATGGACTTGTCGTTGCCACCGAGCATAGCGGTTTCGAAATAGCCGAGATGGGCACTGGCGAATTAACCCACTGGTTTATATGGATCCGATGGGTTCACTTAATAAGCTCGTACCAGGGATGAATAAAGCGTTACGAGAATTATAAACATGGAGTTCCTATTGATTTGAGGTTAATACCGAACGGGAACATTTGTCGATCATGCTTCACATAGAGT'

# Perform approximate string matching between the spike sequence and another set of sequences called ASVs (amplicon sequence variants).
# The result is stored in the variable 'spike_ix'.
spike_ix = agrep(spike, asv)

# Select rows from the 'taxa' variable based on the indices in 'spike_ix'.
# The selected rows represent the spike ASVs.
spiketab = seqtab[spike_ix, ]

# Excluding non-annotated ASVs from the 'taxa' variable, including the spike ASVs.
# The indices of non-annotated ASVs are stored in 'ix_taxa_incomplete'.
ix_taxa_incomplete = which(!complete.cases(taxa[,2]))

# Select rows from 'seqtab' based on the indices of non-annotated ASVs.
# The selected rows represent the incomplete ASVs.
incompletetab = seqtab[ix_taxa_incomplete, ]

# Calculate the column sums of 'incompletetab' and store the result in a data frame called 'incomplete_abundance'.
incomplete_abundance = data.frame(colSums(incompletetab)) %>%
  dplyr::rename("incomplete_abundance" = colSums.incompletetab.)

# Exclude ASVs belonging to the 'Metazoa' group from the 'taxa' variable.
# The indices of ASVs belonging to 'Metazoa' are stored in 'ix_taxa_metazoa'.
ix_taxa_metazoa = which(taxa[,4] == 'Metazoa')

# Select rows from 'seqtab' based on the indices of 'Metazoa' ASVs.
# The selected rows represent the ASVs belonging to 'Metazoa'.
metazoatab = seqtab[ix_taxa_metazoa, ]

# Calculate the column sums of 'metazoatab' and store the result in a data frame called 'metazoa_abundance'.
metazoa_abundance = data.frame(colSums(metazoatab)) %>%
  dplyr::rename("metazoa_abundance" = colSums.metazoatab.)

# Calculate the column sums of 'spiketab' and store the result in a data frame called 'spike_abundance'.
spike_abundance = data.frame(colSums(spiketab)) %>%
  dplyr::rename("spike_abundance" = colSums.spiketab.)

# Calculate the column sums of 'spiketab' and store the result in a data frame called 'read_abundance'.
# Perform additional calculations on the 'read_abundance' data frame.
read_abundance = data.frame(colSums(spiketab)) %>%
  dplyr::rename("spike_abundance" = colSums.spiketab.) %>%
  mutate(tot.readcount = colSums(seqtab)) %>%
  mutate(reads.nospike = tot.readcount - spike_abundance) %>%
  mutate(reads.nospikenometazoa = tot.readcount - spike_abundance - metazoa_abundance$metazoa_abundance) %>%
  mutate(reads.nospikenometazoanoincomplete = tot.readcount - spike_abundance - metazoa_abundance$metazoa_abundance - incomplete_abundance$incomplete_abundance) %>%
  mutate(NGI_sample_ID_18S = row.names(.))

# Update the 'NGI_sample_ID_18S' column in 'spike_abundance' data frame with row names.
spike_abundance$NGI_sample_ID_18S = row.names(spike_abundance)

# Set the row names of the 'metadata' variable to the values in the 'NGI_sample_ID_18S' column.
row.names(metadata) = metadata$NGI_sample_ID_18S

# Join the 'metadata' and 'spike_abundance' data frames based on matching values in 'spike_abundance' column.
metadata = metadata %>%
  left_join(spike_abundance)

# Update the 'meta_ix' variable to contain only the indices where 'spike_abundance' values are greater than 0.
meta_ix = intersect(meta_ix, which(metadata$spike_abundance > 0))

# Select columns from 'seqtab' based on the updated 'meta_ix'.
seqtab = seqtab[, meta_ix]

# Update the 'metadata' variable to contain only the rows corresponding to the updated 'meta_ix'.
metadata = metadata[meta_ix,]

# Calculate the "spike_factor" based on the ratio of 'spike_per_liter_18S' to 'spike_abundance' for each row in 'metadata'.
# Multiply 'seqtab' with the computed "spike_factor" values and store the result in 'seqtab_spike_norm'.
metadata = metadata %>%
  mutate("spike_factor" = spike_per_liter_18S/spike_abundance)

seqtab_spike_norm = seqtab %*% diag(metadata$spike_factor)

# Update the column names of 'seqtab_spike_norm' with the sample IDs from 'metadata'.
colnames(seqtab_spike_norm) = metadata$NGI_sample_ID_18S

# Save 'seqtab_spike_norm' and 'metadata' as tab-separated files.
write.table(seqtab_spike_norm, file = "data/18S/seqtab_18S_spike_normalised.txt", sep = "\t")
write.table(metadata, file = "data/18S/spike_normalized_metadata.txt", sep = "\t")