####################################################################################################################
#
#   This script performs permutation scans by permuting the genoprobs object.
#     Parallel runs are based on the permutation list on line 154
#
#
#
#   Input (Mainly for shell script; see qtl2_manualPerm_one.sh):
#      1: viewer_data:   Path to the qtl viewer .RData 
#      2: dataset_expr:  Which dataset.* and data type to use. 
#                        Example - dataset.islet.proteins|data (if only data is the expression matrix) or dataset.islet.proteins|data|rz (if 'data' is a list and you want to select the rankz expression matrix)
#      3: int_name:     (Optional) A string with interactive variables in samples dataframe separated by '|'. Example - 'sex' or 'sex|batch'. 'NA' if not used
#      4: num_cores:     Number of cores to run
#      5: id:            The column name in expression matrix that you want to run permutations for
#      6: seed:          Seed number for consistent permutations
#      7: n_perm:        Number of permutations to run.
#      8: chunk_number:  (Optional) Numeric value of the chunck number. 'NA' if not used
#      9: chunk_size:    (Optional) Numeric value of the chunck size. Should be consistent. 'NA' if not used
#
#
#
#   Output: 
#       1: Matrix containing LOD scores for each permutation run of the phenotype as .rds file
#
#
#
#
#   Authors: Duy Pham
#   Date:    September 17, 2019
#   E-mails: duy.pham@jax.org
####################################################################################################################

### Load required library packages
options(stringsAsFactors = FALSE)
library(qtl2) 







### Command line arguments / Variables to change. (See 'Input' above)
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
stopifnot(c('genoprobs', 'K', dataset_expr[1]) %in% ls())






### Get expression data
ds   <- get(dataset_expr[1])
expr <- ds[[dataset_expr[2]]]

if(length(dataset_expr) == 3){expr <- expr[[dataset_expr[3]]]}
stopifnot(id %in% colnames(expr))











# Get covariate matrix and make num_cores numeric
covar     <- ds$covar.matrix
num_cores <- as.numeric(num_cores)













### Setting interactive matrix
if(!tolower(int_name) == 'na'){

   # Get samples dataframe 
   samples <- ds$annot.samples
   stopifnot(strsplit(int_name, split = '|', fixed = TRUE)[[1]] %in% colnames(samples))
  
   int_covar <- covar[,grep(int_name, colnames(covar)), drop = FALSE]
}else{
   int_covar <- NULL  
}
print(int_covar)








### For permutation

# Setting seed
set.seed(as.numeric(seed))


# Save genoprobs to another object and store original sample ids
gp <- genoprobs
orig_id <- dimnames(gp[[1]])[[1]]


# Generate permutations for genoprobs
n_perm <- as.numeric(n_perm)
n_samp <- length(orig_id)
n_mark <- sum(unlist(lapply(genoprobs, function(x) length(dimnames(x)[[3]]))))
perms  <- lapply(1:n_perm, function(x) orig_id[sample(x = 1:n_samp, size = n_samp, replace = FALSE)])




# Parallel perms
perm.rng <- 1:length(perms)

if(!tolower(chunk_number) %in% c('na')){
   chunk_number <- as.numeric(chunk_number)
   chunk_size   <- as.numeric(chunk_size)
   
   max_col = length(perms)
   perm.rng = ((chunk_number - 1) * chunk_size + 1):(chunk_number * chunk_size)
  
   if(perm.rng[length(perm.rng)] > max_col){
      perm.rng = perm.rng[1]:max_col
   }
   print(paste("Mapping columns:", perm.rng[1], "to", perm.rng[length(perm.rng)]))
}










### Permutation scans begin
scan1_output <- list()
for(i in 1:length(perm.rng)){
  
    gp <- lapply(genoprobs, function(x){dimnames(x)[[1]] <- perms[[perm.rng[i]]]; x})
    attributes(gp) <- attributes(genoprobs)

    scan1_output[[i]] <- scan1(genoprobs = gp, 
                               pheno     = expr[,id,drop = FALSE],
                               kinship   = K,
                               addcovar  = covar,
                               intcovar  = int_covar,
                               cores     = num_cores)
    
    print(i)
}
scan1_output <- do.call(cbind, scan1_output)
colnames(scan1_output) <- paste0(colnames(scan1_output),'_perm_',perm.rng)














### Save permutation results
if(tolower(int_name) == 'na'){
   if(tolower(chunk_number) == 'na'){
      saveRDS(scan1_output, paste0(id, '_additive_permutations.rds'))
   }else{
      saveRDS(scan1_output, paste0(id, '_additive_perm_chunk_',chunk_number,'.rds'))
   }
  
}else{
  
   if(tolower(chunk_number) == 'na'){
      saveRDS(scan1_output, paste0(id, '_', int_name, '_int_permutations.rds'))
   }else{
      saveRDS(scan1_output, paste0(id, '_', int_name, '_int_perm_chunk_', chunk_number, '.rds'))
   }
}


print('Finished')

