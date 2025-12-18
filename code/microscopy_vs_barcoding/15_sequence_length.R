library(tidyverse)
library(phyloseq)

# Paths
asv_seq_file = 'data/18S/ASVs_18S.txt' ## ASV sequences
phyloseq_file = 'data/18S/phyloseqs.RData' 

## Read preproccesing output
load(phyloseq_file)
asv = read.table(asv_seq_file)$x

tax.tab = tax_table(phyloseq_mb)

asv_id = rownames(tax.tab)

id = as.numeric(str_remove(asv_id, "ASV"))

selected_asv = asv[id]

max(nchar(selected_asv))
min(nchar(selected_asv))
mean(nchar(selected_asv))
