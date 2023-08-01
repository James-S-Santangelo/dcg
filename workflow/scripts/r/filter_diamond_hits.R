library(tidyverse)

protein_blast <- read_delim(snakemake@input[["hits"]], delim = "\t")

protein_blast_filt <- protein_blast %>%
  mutate(product = str_remove(pattern = "^[^\\s]*\\s", stitle)) %>%
  dplyr::select(-stitle) %>%
  group_by(qseqid) %>%
  filter(evalue == min(evalue)) %>%
  filter(pident == max(pident)) %>%
  filter(!grepl("hypothetical", product, ignore.case = TRUE)) %>%
  filter(!grepl("uncharacterized", product, ignore.case = TRUE)) %>%
  filter(!grepl("unnamed", product, ignore.case = TRUE)) %>%
  sample_n(size = 1) %>%
  mutate(note = sprintf("Product pulled from best Diamond hit (pident = %s)", pident)) %>%
  dplyr::select(qseqid, product, note)

write_delim(protein_blast_filt, snakemake@output[["out"]], delim = "\t")
