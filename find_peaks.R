####################################################################################################################
#
#   This script performs the "finding peaks" analysis of qtl2 scan1.
#
#
#
#   Input:
#       1:            Path or filename of qtl viewer .RData 
#       2: scan1:     Path or filename of the output from scan1 (LOD matrix)
#       3: dataset:   Which dataset in the qtl viewer to save lod summary to
#       4: thr:       threshold to find peaks function.
#       5: num_cores: number of cores to run
#       6: type_scan: 'additive' or blank
#       7: type_data: which expression data, 'norm', 'rankz', 'raw' or blank if 6 is blank
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
options(stringsAsFactors = FALSE)
library(qtl2) 
library(dplyr)



### Command line arguments / variables to change
args <- commandArgs(trailingOnly = TRUE)
load(args[1])
scan1_mat <- readRDS(args[2])
dataset   <- get(args[3])
thr       <- as.numeric(args[4])
num_cores <- as.numeric(args[5])
type_scan <- args[6]
type_data <- args[7]

stopifnot("map" %in% ls())
map       <- get("map")








### Perform find_peaks
peaks <- find_peaks(scan1_mat, map = map, threshold = thr, cores = num_cores)










### Formatting for QTL Viewer
peaks <- peaks %>%
               select(-lodindex) %>%                         # Remove lodindex column
               left_join(markers, by = c('chr', 'pos')) %>%  # Get marker ids
               rename(annot.id  = lodcolumn,
                      qtl.chr   = chr,
                      qtl.pos   = pos,
                      marker.id = marker) %>%                # Renaming columns
               select(annot.id, marker.id, lod, qtl.chr, qtl.pos)       # Reorder data frame


print(head(peaks))



### BLUP scan at each QTL
if(type_scan == 'additive'){
    
   stopifnot(c('K','genoprobs') %in% ls())
   peaks = cbind(peaks, matrix(0, nrow = nrow(peaks), ncol = 8,
                               dimnames = list(NULL, LETTERS[1:8])))
   expr  = dataset[[type_data]]
   covar = dataset$covar

   print(head(peaks))
   for(i in 1:nrow(peaks)){
       
       chr     = peaks$qtl.chr[i]
       mkr     = peaks$marker.id[i]
       annot   = peaks$annot.id[i]
       gp      = genoprobs[,chr]
       gp[[1]] = gp[[1]][,,mkr,drop=FALSE]
       
       
       blup    = scan1blup(genoprobs = gp, pheno = expr[,annot, drop = FALSE],
                           kinship = K[[chr]], addcovar = covar, cores = num_cores)
                           
       
       peaks[i, LETTERS[1:8]] <- blup[1,1:8]
      
    }
    
    
}



dataset$lod.peaks[[type_scan]] <- peaks
assign(args[3], dataset)


save(list = ls()[ls() %in% c(grep('dataset.', ls(), value = TRUE),
                                  'K',
                                  'genoprobs',
                                  'map',
                                  'markers',
                                  'ensemble.version')], file = args[1])





