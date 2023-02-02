# Script to generate pairwise haplotype dotplots

###############
#### SETUP ####
###############

# Load required packages
library(tidyverse)
library(pafr)

# Source raw code for pafr dotplot() function so it can be modified below
source('https://raw.githubusercontent.com/dwinter/pafr/master/R/dotplot.R')

# Minor tweaks to dotplot() function for nicer plots
dotplot_mod <- function(ali, order_by = c("size", "qstart", "provided"),
                    label_seqs = FALSE, dashes = TRUE, ordering = list(),
                    alignment_colour = "black", xlab = "query", ylab = "target",
                    line_size=2) {
    by <- match.arg(order_by)
    if (by == "provided") {
        check_ordering(ali, ordering)
        ali <- ali[ali$qname %in% ordering[[1]] & ali$tname %in% ordering[[2]],]
    }
    seq_maps <- order_seqs(ali, by, ordering)
    ali <- add_pos_in_concatentaed_genome(ali, seq_maps)
    #now build the plot, first the aligments, easiest to do this +ve strand
    #first, then -ve
    p <- ggplot() +
      geom_segment(data = ali[ali$strand=="+",], 
         aes_string(x = "concat_qstart", xend = "concat_qend", y = "concat_tstart", yend = "concat_tend"),
         size=line_size, colour=alignment_colour) +
      geom_segment(data = ali[ali$strand=="-",],
        aes_string(x = "concat_qend", xend = "concat_qstart", y = "concat_tstart", yend = "concat_tend"),
        size = line_size, colour = alignment_colour) +
      coord_equal() +
      scale_x_continuous(xlab, labels = Mb_lab) + 
      scale_y_continuous(ylab, labels = Mb_lab)
  if (dashes) {
      p <- p + geom_hline(yintercept = c(seq_maps[["tmap"]], sum(unique(ali$tlen))), linetype = 3) +
               geom_vline(xintercept = c(seq_maps[["qmap"]], sum(unique(ali$qlen))), linetype = 3)
  }
  if (label_seqs) {
    qname_df <- dotplot_name_df(seq_maps[["qmap"]], seq_maps[["qsum"]])
    tname_df <- dotplot_name_df(seq_maps[["tmap"]], seq_maps[["tsum"]])
    p <- p + geom_text(data = qname_df, aes_string(label = "seq_name", x = "centre", y = "0"),
                       vjust = 1, angle = 90) +
             geom_text(data = tname_df,
                       aes_string(label = "seq_name", x = "0", y = "centre"), vjust = 0)
  }
  # We want to be able to annotated the dotpot with data in BED format. Adding
  # arugments to this fxn wouldbe pretty unwieldy, so we want to take advantage
  # of ggplots overloaded "+" to add the annotatoins. To oders the annotations
  # in that same way as the dotplot object, we need this object to include a
  # functoin for mapping chromosome positoins to concatenated genome positoins
  # in the dotlot
  p$seq_map_fxn <- function(bed, query=TRUE, ...){
      map_n <-  if (query) 1 else 3
      seq_map <- seq_maps[[map_n]]
      check_chroms <-  bed[["chrom"]] %in% names(seq_map)
      if( !(all(check_chroms)) ){
          if( !(any(check_chroms))) {
              stop("None of the chromosomes represented this bed file are part of the dotplot")
          } else {
              bed <- bed[bed$chrom %in% names(seq_map),]
              missing <- unique( bed[["chrom"]] [!check_chroms])
              msg <- paste(length(missing), 
                           "of the chromosomes in this bed file are not part of the dotplot:\n  ", 
                           paste(missing, collapse=", "))
              warning(msg, call.=FALSE)
          }
      }
      data.frame(istart=bed[["start"]] + seq_map[ bed[["chrom"]] ], 
                 iend = bed[["end"]] + seq_map[ bed[["chrom"]] ],
                len = seq_maps[[map_n + 1]] 
      )
  }
  p
}

# Function to load PAF file. Leverages `pafr` package
load_filter_paf <- function(path, target_hap, query_hap, mapq = 0, alen = 0){
    
    paf <- pafr::read_paf(path) %>%
        dplyr::select(qname:mapq)

    scaff_regex <- '^Scaffold_([1-9]|10)(?=_)'
    extract <- '(?<=Scaffold_)([1-9]|10)(?=_)'
    
    paf_mod <- paf %>%
        filter(mapq >= mapq & alen >= alen) %>%
        mutate(pid = nmatch / alen) %>%
        mutate(queryhap = query_hap,
               targethap = target_hap) %>%
    
        # Get only mappings with target and query on first 10 scaffolds
        filter(str_detect(qname, scaff_regex) & str_detect(tname, scaff_regex)) %>%
    
        # Rename target and query sequences with haplotype info
        mutate(qname = paste0(queryhap, '_S', str_extract(qname, extract)),
               tname = paste0(targethap, '_S', str_extract(tname, extract)))
    
    return(paf_mod)
}

# Function to plot dotplot for particular haplotypes
plot_dotplot <- function(df, target_hap, query_hap){
    
    plt <- df %>%
        filter(targethap == target_hap & queryhap == query_hap) %>%
        dotplot_mod(., order_by = "qstart", label_seqs=TRUE) + 
        ggtitle(sprintf("Target = %s, Query = %s", target_hap, query_hap)) + 
        theme_classic() +
        theme(axis.title = element_text(size = 17),
              axis.text = element_text(size = 14),
              plot.title = element_text(size = 20, face = 'bold'))
    
    ggsave(plt, filename = dotplot_out, device = 'pdf', width = 10, height = 10, 
            units = 'in', dpi = 300)
}

##########################
#### GENERATE DOTPLOT ####
##########################

# Define inputs and outputs
paf_in <- snakemake@input[[1]]
dotplot_out <- snakemake@output[[1]]

# Get haplotype names from input filepath
split <- str_split(basename(tools::file_path_sans_ext(paf_in)), '-', simplify = T)
target <- split[1]
query <- split[2]

# Load PAF file
paf <- load_filter_paf(paf_in, target_hap = target, query_hap = query)

plot_dotplot(paf, target_hap = target, query_hap = query)





