library(GENESPACE)

###############################################
# -- change paths to those valid on your system
genomeRepo <- snakemake@input['dir'] 
wd <- snakemake@output[[1]]
path2mcscanx <- "/opt/MCScanX-master"
###############################################

sink(snakemake@log[[1]])
# -- parse the annotations to fastas with headers that match a gene bed file
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = c("Pall", "Occ"),
  presets = "none",
  gffIdColum = "ID",
  headerEntryIndex = 1,
  headerSep = " ",
  gffStripText = "",
  headerStripText = "",
  genespaceWd = wd,
)

# -- initalize the run and QC the inputs
gpar <- init_genespace(
  wd = wd,
  useHOGs = TRUE,
  orthofinderInBlk = TRUE,
  path2mcscanx = path2mcscanx,
  nCores = snakemake@threads)

# -- accomplish the run
out <- run_genespace(gpar)
sink()
