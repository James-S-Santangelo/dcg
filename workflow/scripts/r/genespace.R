library(GENESPACE)

###############################################
# -- change paths to those valid on your system
genomeRepo <- snakemake@input['dir'] 
wd <- snakemake@output[[1]]
path2mcscanx <- "/opt/MCScanX-master"
###############################################

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
  troubleShoot = TRUE,
  genespaceWd = wd,
)

# -- initalize the run and QC the inputs
gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx,
  useHOGs = TRUE,
  orthofinderInBlk = TRUE,
  nCores = snakemake@threads)

# -- accomplish the run
out <- run_genespace(gpar)
