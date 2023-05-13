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
                       vjust = 1, angle = 90, size = 5) +
             geom_text(data = tname_df,
                       aes_string(label = "seq_name", x = "0", y = "centre"), vjust = 0, size = 5)
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

# Function to plot dotplot for particular haplotypes
plot_dotplot <- function(df, ver){
    title <- sprintf("Target = Hap2, Query = Hap1, Version: %s", ver)
    plt <- df %>%
        dotplot_mod(., order_by = "qstart", label_seqs=TRUE) + 
        ggtitle(title) + 
        theme_classic() +
        theme(axis.title = element_text(size = 17),
              axis.text = element_text(size = 17),
              plot.title = element_text(size = 20, face = 'bold'))
    
    ggsave(plt, filename = dotplot_out, device = 'pdf', width = 14, height = 14, 
            units = 'in', dpi = 300)
}

##########################
#### GENERATE DOTPLOT ####
##########################

# Define inputs and outputs
paf_in <- snakemake@input[[1]]
dotplot_out <- snakemake@output[[1]]
ver <- snakemake@wildcards[['ver']]

# Load PAF file
paf <- pafr::read_paf(paf_in) %>% dplyr::select(qname:mapq)

plot_dotplot(paf, ver)





