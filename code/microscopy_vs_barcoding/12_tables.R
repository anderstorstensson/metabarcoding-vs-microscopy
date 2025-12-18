library(tidyverse)
library(phyloseq)
library(cowplot)
library(MicrobiotaProcess)

# File paths 
phyloseq_file = 'data/18S/phyloseqs.RData' 
station_file <- 'data/everything_from_shark/stations_to_extract.txt'

# Read files
stations <- read.table(station_file, sep = '\t', header = TRUE) %>%
  dplyr::select(-Phytoplankton_analysis,
                -nr_of_samples,
                -sea_basin) %>%
  mutate(SHARK = iconv(SHARK, from = "latin1", to = "UTF-8"),
         station_plots = iconv(station_plots, from = "latin1", to = "UTF-8"))

# Load phyloseq objects
load(phyloseq_file)

metadata = data.frame(sample_data(phyloseq_mb)) %>%
  mutate(station_name = gsub("R\u0081NE\u0081-2", "RÅNEÅ-2", station_name)) %>%
  mutate(station_name = gsub("SL€GG…", "SLÄGGÖ", station_name)) %>%
  mutate(station_name = gsub("\u008117", "Å17", station_name)) %>%
  mutate(station_name = gsub("R\u0081NE\u0081-1", "RÅNEÅ-1", station_name)) %>%
  mutate(station_name = gsub("SR3", "SR3/C24", station_name)) %>%
  left_join(stations, by = c("station_name" = "SHARK")) %>%
  mutate(percent_phytoplankton = (colSums(otu_table(phyloseq_mb))/sample_data(phyloseq_mb)$reads.nospikenometazoa)*100)

station_table = metadata %>%
  group_by(station_plots) %>%
  summarise(station_name = unique(station_name),
            sea_basin = unique(sea_basin),
            n_sample = length(unique(NGI_sample_ID_18S)),
            salinity = round(mean(salinity), 1),
            dna_conc = round(mean(sample_DNA_concentration), 1),
            time = paste(month.abb[min(sampling_month)], month.abb[max(sampling_month)], sep="-"),
            n_month = length(unique(sampling_month)),
            percent_phytoplankton = round(mean(percent_phytoplankton), 0)) %>%
  arrange(sea_basin)

print(paste("Total amount of sequences:", sum(colSums(otu_table(phyloseq_mb)))))
print(paste("Max sequences per sample:", max(colSums(otu_table(phyloseq_mb)))))
print(paste("Min sequences per sample:", min(colSums(otu_table(phyloseq_mb)))))
print(paste("Number of ASVs:", nrow(phyloseq::tax_table(phyloseq_mb))))

print(paste("Number of microscopy taxa:", nrow(phyloseq::tax_table(phyloseq_abutab))))
print(paste("Max cells per sample:", max(colSums(otu_table(phyloseq_abutab)))))
print(paste("Min cells per sample:", min(colSums(otu_table(phyloseq_abutab)))))

write.table(station_table, "data/tables/station_summary.txt", sep = "\t", row.names = FALSE)
