####################################################################################################################
#
#   This script gathers all qtl2 scan1 matrices (chunks) and cbinds them together.
#
#
#
#
#   Input:
#       1: pattern:   Pattern of the qtl2 chunk file name
#       2: out_file:  File name to save the concatenated chunks without the '.rds'
#
#
#   Output: 
#       1: Matrix of all chunks cbind together
#
#
#
#
#   Authors: Duy Pham, Andrew Deighan, & Isabela Gyuricza
#   Date:    October 31, 2018
#   E-mails: duy.pham@jax.org, andrew.deighan@jax.org, & isabela.gyuricza@jax.org
####################################################################################################################


### Load required library packages
options(stringsAsFactors = FALSE)
library(data.table)




### Command line arguments / variables to change
# 1.) pattern: pattern without the _chunk_#.rds
# 2.) chunk start: chunk start
# 3.) chunk_end: chunk end
# 4.) out_file: name to save the combined chunk file
args <- commandArgs(trailingOnly = TRUE)
pattern     <- args[1]
chunk_start <- as.numeric(args[2])
chunk_end   <- as.numeric(args[3])
out_file    <- args[4]





### Read in all qtl2 chunk file name and cbind them together
temp <- list()
for(i in chunk_start:chunk_end){
    temp[[i]] <- readRDS(paste0(pattern,'_chunk_',i,'.rds'))
}

temp <- do.call(cbind, temp)




### Save the matrix as .rds
saveRDS(temp, file = paste0(out_file,'.rds'))
