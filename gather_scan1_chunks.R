####################################################################################################################
#
#   This script gathers all qtl2 scan1 matrices (chunks) and combines them together.
#
#
#
#
#   Input:
#       1: pattern:     Pattern of the qtl2 chunk file name including 'chunk'. Ex. islet_proteins_additive_scan_chunk'
#       2: chunk_start: Starting chunk number
#       3: chunk_end:   Ending chunk number
#       4: func:        Function used to bind all chunks together. Ex. 'rbind' or 'cbind'
#       2: out_file:    File name to save the concatenated chunks without the '.rds'
#
#
#   Output: 
#       1: Matrix of all chunks combined
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
args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)==0){
    print("No arguments supplied.")
}else{
    for(i in 1:length(args)){
        a <- strsplit(args[i],split = '=', fixed = TRUE)[[1]]
        assign(a[1],a[2])
    }
}





### Convert the chunk starting and ending numbers to numeric
chunk_start = as.numeric(chunk_start)
chunk_end = as.numeric(chunk_end)







### Read in all chunk files and combine them together
temp <- list()
for(i in chunk_start:chunk_end){
    	temp[[i]] <- readRDS(paste0(pattern,'_chunk_',i,'.rds'))
}

temp <- do.call(func, temp)








### Save the matrix as .rds
saveRDS(temp, file = paste0(out_file,'.rds'))
