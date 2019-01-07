####################################################################################################################
#
#   This script performs qtl2 additive scan
#
#   Notes*:
#         The QTL scan can be run in 'chunks' or all at once.
#             
#             Example for chunk size/number: 
#             Suppose there are 5433 phenotype columns. If chunk size is 1000, then there should be 6 different chunk_number,
#                  with the chunk_number value being 1,2,3,4,5, or 6 to get the columns:
#                       1-1000,1001-2000,2001-3000,3001-4000,4001-5000,5001-5433, respectively.
#
#         If you do not want to run in chunks, leave chunk_number and chunk_size blank
#
#
#
#
#   Input:
#       1: Input file:    Path to the qtl viewer .RData 
#       2: dataset:       Which dataset.* to use
#       3: expression:    Which transformation of the phenotype to use. "norm" or "rankz"
#       4: int_term:      a string with interactive variables seperated by ','. Ex. 'sex' or 'sex,batch'
#       5: num_cores:     Number of cores to run
#       6: chunk_number:  Numeric value of the chunk number.  Can be left blank
#       7: chunk_size:    Numeric value of chunk size. Should be consistent.  Can be left blank
#
#   Output: 
#       1: Matrix containing LOD scoress for each of the phenotype that was given to scan1 at each marker.
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
library(qtl2) 







### Command line arguments / Variables to change
# 1: input.file:    Path to the qtl viewer .RData 
# 2: dataset:       Which dataset.* to use
# 3: expression:    Which transformation of the phenotype to use. "norm" or "rankz"
# 4: int_term:      a string with interactive variables in samples dataframe separated by '|'.     Ex. 'sex' or 'sex|batch'
# 5: num_cores:     Number of cores to run
# 6: chunk_number:  Numeric value of the chunk number. Can be left blank
# 7: chunk_size:    Numeric value of chunk size. Should be consistent. Can be left blank
args = commandArgs(trailingOnly = TRUE)

load(args[1])
dataset       <- get(args[2])
expr          <- dataset[[args[3]]]
int_term      <- args[4]
num_cores     <- as.numeric(args[5])
chunk_number  <- as.numeric(args[6])
chunk_size    <- as.numeric(args[7])


# Get covariate and samples matrix
covar   <- dataset$covar
samples <- dataset$samples







### Check to see if required data are loaded in global environment
stopifnot(strsplit(int_term, split = '|', fixed = TRUE)[[1]] %in% colnames(samples))
stopifnot(c("genoprobs", "K", args[2]) %in% ls())






### Get interactive matrix
int_covar <- covar[,grep(int_term, colnames(covar)), drop = FALSE]






### Run qtl scan on a subset on a subset of the expr matrix if chunk_size and chunk_number is given
pheno.rng <- 1:ncol(expr)

if(!is.na(chunk_number)){
   max_col = ncol(expr)
   pheno.rng = ((chunk_number - 1) * chunk_size + 1):(chunk_number * chunk_size)
 
   if(pheno.rng[length(pheno.rng)] > max_col){
    
      pheno.rng = pheno.rng[1]:max_col
 
   }
  
   print(paste("Mapping columns:", pheno.rng[1], "to", pheno.rng[length(pheno.rng)]))

}












### qtl2 scan1 
scan1_output <- scan1(genoprobs = genoprobs, 
                      pheno     = expr[,pheno.rng],
                      kinship   = K,
                      addcovar  = covar,
                      intcovar  = int_covar,
                      cores     = num_cores)











### Save scan1 output
output_file <- sub('dataset.', '', args[2], fix = TRUE)
output_file <- sub('.', '_', output_file, fix = TRUE)
output_int  <- gsub('|','_',int_term)

if(is.na(chunk_number)){
   saveRDS(scan1_output, paste0(output_file, '_', output_int, '_int_scan.rds'))
}else{
   saveRDS(scan1_output, paste0(output_file, '_', output_int, '_int_scan_chunk_', chunk_number,'.rds'))
}
