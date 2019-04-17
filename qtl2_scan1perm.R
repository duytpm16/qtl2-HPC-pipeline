####################################################################################################################
#
#   This script performs qtl2 permutation scan
#
#   Notes*:
#         The QTL scan can be run in 'chunks' or all at once.
#             
#             Example for chunk size/number: 
#             Suppose there are 5433 phenotype columns. If chunk size is 1000, then there should be 6 different chunk_number,
#                  with the chunk_number value being 1,2,3,4,5, or 6 to get the columns:
#                       1-1000,1001-2000,2001-3000,3001-4000,4001-5000,5001-5433, respectively.
#
#         If you do not want to run in chunks, make chunk_number and chunk_size 'NA' or 'na'
#
#
#
#
#   Input:
#      1: input.file:    Path to the qtl viewer .RData 
#      2: dataset:       Which dataset.* to use
#      3: int_name:      (Optional) A string with interactive variables in samples dataframe separated by '|'. Ex. 'sex' or 'sex|batch'. 'NA' or 'na' if not used
#      4: num_cores:     Number of cores to run
#      5: chunk_number:  (Optional) Numeric value of the chunk number. 'NA' or 'na' if not used
#      6: chunk_size:    (Optional) Numeric value of chunk size. Should be consistent.  'NA' or 'na' if not used
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





### Load viewer data
load(viewer_data)






### Check to see if required data are loaded in global environment
dataset_expr <- strsplit(dataset_expr, '|',fixed = TRUE)[[1]]
print(dataset_expr)
stopifnot(c('genoprobs', 'K', dataset_expr[1]) %in% ls())

ds    <- get(dataset_expr[1])
expr  <- ds[[dataset_expr[2]]]

if(length(dataset_expr) == 3){
   expr <- expr[[dataset_expr[3]]]
}












int_term  <- int_name
num_cores <- as.numeric(num_cores)
perm_run  <- as.numeric(perm_run)



# Get covariate matrix
covar   <- ds$covar


















### Get interactive matrix
if(!int_term %in% c('NA','na')){
 
   # Get samples dataframe
   samples <- ds$samples
   stopifnot(strsplit(int_term, split = '|', fixed = TRUE)[[1]] %in% colnames(samples))

   int_covar <- covar[,grep(int_term, colnames(covar)), drop = FALSE]
   print(int_covar)
}else{

   int_covar <- NULL

}









### Run qtl scan on a subset on a subset of the expr matrix if chunk_size and chunk_number is given
pheno.rng <- 1:ncol(expr)

if(!chunk_number %in% c('NA','na')){
   chunk_number <- as.numeric(chunk_number)
   chunk_size   <- as.numeric(chunk_size)

   max_col = ncol(expr)
   pheno.rng = ((chunk_number - 1) * chunk_size + 1):(chunk_number * chunk_size)

   if(pheno.rng[length(pheno.rng)] > max_col){

      pheno.rng = pheno.rng[1]:max_col

   }

   print(paste("Mapping columns:", pheno.rng[1], "to", pheno.rng[length(pheno.rng)]))

}












### qtl2 scan1 
scan1_output <- scan1perm(genoprobs = genoprobs, 
                      	  pheno     = expr[,pheno.rng,drop = FALSE],
                          kinship   = K,
                          addcovar  = covar,
			  intcovar  = int_covar,
                          cores     = num_cores,
                          n_perm    = perm_run)













### Save scan1 output
output_file <- sub('dataset.', '', dataset_expr[1], fixed = TRUE)
output_file <- sub('.', '_', output_file, fixed = TRUE)
output_int  <- gsub('|','_',int_term, fixed = TRUE)


if(int_term == 'NA'){
   if(chunk_number == 'NA'){
      saveRDS(scan1_output, paste0(output_file, '_perm_additive_scan.rds'))
   }else{
      saveRDS(scan1_output, paste0(output_file, '_perm_additive_scan_chunk_',chunk_number,'.rds'))
   }

}else{

   if(chunk_number == 'NA'){
      saveRDS(scan1_output, paste0(output_file, '_', output_int, '_perm_int_scan.rds'))
   }else{
      saveRDS(scan1_output, paste0(output_file, '_', output_int, '_perm_int_scan_chunk_', chunk_number,'.rds'))
   }
}

