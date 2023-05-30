library(tidyverse)

# Function to load results of BLASTING linkage markers against the haplotypes
load_marker_blast <- function(path){
 
  path_split <- str_split(basename(path), pattern = "_", simplify = TRUE)
  marker_pop <- path_split[2][1]
  ref <- path_split[1][1]

  # Column names
  cols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "qlen", "sstart", "send", "slen", "evalue",
            "bitscore", "qcovs", "qcovhsp")

  # Read in BLAST results
  markers_blast <- read.delim(path, sep = "\t", col.names = cols) %>%
    mutate(qseqid = as.character(qseqid)) %>%
    mutate(pop = marker_pop,
           ref = ref) %>%
    mutate(sseqid = case_when(ref == "TrRv5" ~ str_extract(sseqid, "(?<=\\|).*(?=\\|)"),
                              ref == "utm" ~ sseqid)) %>%

  return(markers_blast)
}

# Inputs
markers_blast_paths <- snakemake@input[1:4]
DG_names_path <- snakemake@input[[5]]
SG_names_path <- snakemake@input[[6]]
DG_genMap_path <- snakemake@input[[7]]
SG_genMap_path <- snakemake@input[[8]]

# Outputs
gm_kar_DG <- snakemake@output[[1]]
match_gm_links_DG <- snakemake@output[[2]]
nomatch_gm_links_DG <- snakemake@output[[3]]
markerPos_DG <- snakemake@output[[4]]
gm_kar_SG <- snakemake@output[[5]]
match_gm_links_SG <- snakemake@output[[6]]
nomatch_gm_links_SG <- snakemake@output[[7]]
markerPos_SG <- snakemake@output[[8]]

# Load results from BLAST searches
genMap_df <- purrr::map_dfr(markers_blast_paths, load_marker_blast)

# Load dataframes with names of markers.
DG_names <- read.csv(DG_names_path, col.names = c("qseqid", "name")) %>% mutate(pop = "DG")
SG_names <- read.csv(SG_names_path, col.names = c("qseqid", "name")) %>% mutate(pop = "SG")
all_names <- bind_rows(SG_names, DG_names) %>% mutate(qseqid = as.character(qseqid))

# Load dataframes with genetic positions of markers.
DG_genMap <- read.csv(DG_genMap_path, col.names = c("name", "LG", "cM")) %>% mutate(pop = "DG")
SG_genMap <- read.csv(SG_genMap_path, col.names = c("name", "LG", "cM")) %>% mutate(pop = "SG")
all_genMap <- bind_rows(SG_genMap, DG_genMap)

# Join all dataframes
genMap_df_mod <- genMap_df %>%
  left_join(., all_names, by = c("qseqid","pop")) %>%
  left_join(., all_genMap, by = c("name", "pop"))

# Filter
genome_size <- 968215888
DG_mapLength <- 5057.6
SG_mapLength <- 5815.6

filtered_markers <- genMap_df_mod %>%
  group_by(name, ref) %>%
 
  # Remove cases where the sanem marker maps identically to multiple pos on same chrom
  distinct(qseqid, sseqid, pident, evalue, bitscore, .keep_all = T) %>%
 
  # For each marker and reference genome, get only long alignments with high % ident.
  # Filter for best alignment
  filter(pident >= 95 & length > 175) %>%
  filter(evalue == min(evalue)) %>%
  ungroup() %>%

  # Calculate mid positions of markers on LG and chromosomes
  mutate(lg_pos = round(case_when(pop == "DG" ~ cM / (DG_mapLength / genome_size),
                                  pop == "SG" ~ cM / (SG_mapLength / genome_size))),
         chrom_pos = case_when(sstart < send ~ sstart + (qlen / 2),
                               sstart > send ~ sstart - (qlen / 2))) %>%

  dplyr::select(name, sseqid, pident, length, pop, ref, LG, cM, lg_pos, chrom_pos) %>%

  # Filter if no genetic position or identical duplicate alignments to diff. chromosomes
  filter(!is.na(cM)) %>%
  group_by(ref, name) %>%
  filter(n() == 1) %>%
  ungroup()

###############################
#### DG MAPPING POPULATION ####
###############################

filtered_markers_DG <- filtered_markers %>%
  filter(pop == "DG")

# Create karyotype file for linkage map
LG_karyotype_DG <- filtered_markers_DG %>%
  group_by(LG) %>%
  mutate(end = as.integer(max(lg_pos))) %>%
  ungroup() %>%
  mutate(LG = paste0("LG_", LG)) %>%
  mutate(chr = "chr", x = "-", y = 0, col = "white", label = LG) %>%
  dplyr::select(chr, x, LG, label, y, end, col) %>%
  distinct()
write_delim(LG_karyotype_DG, gm_kar_DG, delim = "\t", col_names = F)

# Keep only markers with alignments to both references
markers_toKeep_DG <- dplyr::intersect(
  filtered_markers_DG %>% filter(ref == "TrRv5") %>% dplyr::select(name),
  filtered_markers_DG %>% filter(ref == "utm") %>% dplyr::select(name)
) %>%
  pull(name)

# Dataframe LG to Chromosome maps
utm_LG_map <- data.frame(
  sseqid = c(
    "Chr01_Occ", "Chr01_Pall", "Chr02_Occ", "Chr02_Pall", "Chr03_Occ", "Chr03_Pall",
    "Chr04_Occ", "Chr04_Pall", "Chr05_Occ", "Chr05_Pall", "Chr06_Occ", "Chr06_Pall",
    "Chr07_Occ", "Chr07_Pall", "Chr08_Occ", "Chr08_Pall"
    ),
  corr_LG = c("LG_1", "LG_9", "LG_2", "LG_10", "LG_3", "LG_11", "LG_4", "LG_12",
              "LG_5", "LG_13", "LG_6", "LG_14", "LG_7", "LG_15", "LG_8", "LG_16")
)
TrRv5_LG_map <- data.frame(
  sseqid = c(
    "CM019101.1", "CM019102.1", "CM019103.1", "CM019104.1", "CM019105.1", "CM019106.1",
    "CM019107.1", "CM019108.1", "CM019109.1", "CM019110.1", "CM019111.1", "CM019112.1",
    "CM019113.1", "CM019114.1", "CM019115.1", "CM019116.1"
  ),
  corr_LG = c("LG_1", "LG_9", "LG_2", "LG_10", "LG_3", "LG_11", "LG_4", "LG_12",
              "LG_5", "LG_13", "LG_6", "LG_14", "LG_7", "LG_15", "LG_8", "LG_16")
)

# Create file with links for UTM reference to LG
utm_genMap_links_DG <- filtered_markers_DG %>%
  filter(ref == "utm" & name %in% markers_toKeep_DG) %>%
  mutate(LG = paste0("LG_", LG),
         cspos = chrom_pos, cepos = chrom_pos,
         lspos = lg_pos, lepos = lg_pos) %>%
  left_join(., utm_LG_map, by = "sseqid") %>%
  mutate(color = case_when(
    LG == corr_LG ~ "grey",
    TRUE ~ "red"
  ))

# Create file with links for Griffiths reference to LG
TrRv5_genMap_links_DG <- filtered_markers_DG %>%
  filter(ref == "TrRv5" & name %in% markers_toKeep_DG) %>%
  mutate(LG = paste0("LG_", LG),
         cspos = chrom_pos, cepos = chrom_pos,
         lspos = lg_pos, lepos = lg_pos) %>%
  left_join(., TrRv5_LG_map, by = "sseqid") %>%
  mutate(color = case_when(
    LG == corr_LG ~ "grey",
    TRUE ~ "red"
  ))

# Link files for markers that match and or don"t match correctLG
allLinks_DG <- bind_rows(utm_genMap_links_DG, TrRv5_genMap_links_DG)
matchLinks_DG <- allLinks_DG %>%
  filter(LG == corr_LG) %>%
  dplyr::select(sseqid, cspos, cepos, LG, lspos, lepos)
write_delim(matchLinks_DG, match_gm_links_DG, col_names = FALSE, delim = "\t")

noMatchLinks_DG <- allLinks_DG %>%
  filter(LG != corr_LG) %>%
  dplyr::select(sseqid, cspos, cepos, LG, lspos, lepos)
write_delim(noMatchLinks_DG, nomatch_gm_links_DG, col_names = FALSE, delim = "\t")
 
# Dataframe with marker positions
marker_pos_DG <- utm_genMap_links_DG %>%
  dplyr::select(LG, lspos, lepos)
write_delim(marker_pos_DG, markerPos_DG, col_names = FALSE, delim = "\t")

allLinks_DG %>%
  mutate(is_match = ifelse(color == 'grey', 1, 0)) %>%
  group_by(ref) %>%
  summarise(n = n(),
            matches = sum(is_match),
            prop = matches / n)
  

###############################
#### SG MAPPING POPULATION ####
###############################

filtered_markers_SG <- filtered_markers %>%
  filter(pop == "SG")

markers_toKeep_SG <- dplyr::intersect(
  filtered_markers_SG %>% filter(ref == "TrRv5") %>% dplyr::select(name),
  filtered_markers_SG %>% filter(ref == "utm") %>% dplyr::select(name)
) %>%
  pull(name)

utm_genMap_links_SG <- filtered_markers_SG %>%
  filter(ref == "utm" & name %in% markers_toKeep_SG) %>%
  mutate(LG = paste0("LG_", LG),
         cspos = chrom_pos, cepos = chrom_pos,
         lspos = lg_pos, lepos = lg_pos) %>%
  left_join(., utm_LG_map, by = "sseqid") %>%
  mutate(color = case_when(
    LG == corr_LG ~ "grey",
    TRUE ~ "red"
  ))

# Create file with links for Griffiths reference to LG
TrRv5_genMap_links_SG <- filtered_markers_SG %>%
  filter(ref == "TrRv5" & name %in% markers_toKeep_SG) %>%
  mutate(LG = paste0("LG_", LG),
         cspos = chrom_pos, cepos = chrom_pos,
         lspos = lg_pos, lepos = lg_pos) %>%
  left_join(., TrRv5_LG_map, by = "sseqid") %>%
  mutate(color = case_when(
    LG == corr_LG ~ "grey",
    TRUE ~ "red"
  ))

  # Create karyotype file for linkage map
LG_karyotype_SG <- filtered_markers_SG %>%
  group_by(LG) %>%
  mutate(end = as.integer(max(lg_pos))) %>%
  ungroup() %>%
  mutate(LG = paste0("LG_", LG)) %>%
  mutate(chr = "chr", x = "-", y = 0, col = "white", label = LG) %>%
  dplyr::select(chr, x, LG, label, y, end, col) %>%
  distinct()
write_delim(LG_karyotype_SG, gm_kar_SG, delim = "\t", col_names = F)

# Link files for markers that match and or don"t match correctLG
allLinks_SG <- bind_rows(utm_genMap_links_SG, TrRv5_genMap_links_SG)
matchLinks_SG <- allLinks_SG %>%
  filter(LG == corr_LG) %>%
  dplyr::select(sseqid, cspos, cepos, LG, lspos, lepos)
write_delim(matchLinks_SG, match_gm_links_SG, col_names = FALSE, delim = "\t")

noMatchLinks_SG <- allLinks_SG %>%
  filter(LG != corr_LG) %>%
  dplyr::select(sseqid, cspos, cepos, LG, lspos, lepos)
write_delim(noMatchLinks_SG, nomatch_gm_links_SG, col_names = FALSE, delim = "\t")
 
# Dataframe with marker positions
marker_pos_SG <- utm_genMap_links_SG %>%
  dplyr::select(LG, lspos, lepos)
write_delim(marker_pos_SG, markerPos_SG, col_names = FALSE, delim = "\t")

allLinks_SG %>%
  mutate(is_match = ifelse(color == 'grey', 1, 0)) %>%
  group_by(ref) %>%
  summarise(n = n(),
            matches = sum(is_match),
            prop = matches / n)
