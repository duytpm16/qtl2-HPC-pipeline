####################################################################################################################
#
#   This script performs the "finding peaks" analysis of qtl2 scan1.
#
#
#
#   Input:
#       1:            Path or filename of qtl viewer .RData 
#       2: scan1:     Path or filename of the output from scan1 (LOD table).
#       3: thr:       threshold to find peaks function.
#       4: out_file:  File name (phenotype) to save the output from find peaks function.
#
#
#   Output: 
#       1: Output from find peaks function - QTL peaks above threshold and formatted for QTL Viewer
#
#
#
#
#   Authors: Duy Pham, Andrew Deighan, & Isabela Gyuricza
#   Date:    November 6, 2018
#   E-mails: duy.pham@jax.org, andrew.deighan@jax.org, & isabela.gyuricza@jax.org
#
####################################################################################################################

### Load required library packages
library(qtl2) 
library(dplyr)



### Command line arguments / variables to change
args <- commandArgs(trailingOnly = TRUE)

load(args[1])
stopifnot("map" %in% ls())

scan1_mat <- readRDS(args[2])
map       <- get("map")
thr       <- as.numeric(args[3])
num_cores <- as.numeric(args[4])
type_scan <- args[5]









### Perform find_peaks
peaks <- find_peaks(scan1_mat, map = map, threshold = thr, cores = num_cores)










### Formatting for QTL Viewer
peaks <- peaks %>%
               select(-lodindex) %>%                         # Remove lodindex column
               left_join(markers, by = c('chr', 'pos')) %>%  # Get marker ids
               rename(annot.id  = lodcolumn,
                      qtl.chr   = chr,
                      qtl.pos   = pos
                      marker.id = marker) %>%                # Renaming columns
               select(annot.id, marker.id, lod, qtl.chr, qtl.pos)       # Reorder data frame






### BLUP scan at each QTL
if(args[5] == 'additive'){
    
   stopifnot('genoprobs' %in% ls())
   lod.peaks = cbind(peaks, matrix(0, nrow = nrow(lod.peaks), ncol = 8,
                                   dimnames = list(NULL, LETTERS[1:8])))
   
   for(i in 1:nrow(peaks)){
       
       chr  = peaks$qtl.chr[i]
       mkr  = peaks$marker.id[i]
       gene = peaks$annot.id[i]
       gp      = genoprobs[,chr]
       gp[[1]] = gp[[1]][,,mkr,drop=FALSE]
       
       
       blup    = scan1blup(genoprobs = gp, pheno = expr[,gene, drop = FALSE],
                           kinship = K[[chr]], addcovar = covar, cores = num_cores)
                           
       
       peaks[i, LETTERS[1:8]] <- blup[1,1:8]
      
    }
    
    
}




dataset$lod.peaks[[type_scan]] <- peaks
assign(args[2], dataset)


rm(list = ls()[!(ls() %in% c(grep('dataset.', ls(), value = TRUE),
                                  'K',
                                  'genoprobs',
                                  'map',
                                  'markers',
                                  'ensemble.version'))])


save.image(file = args[1])



