####################################################################################################################
#
#   This script is used to generate SNP association plots
#
#
#
#   Input:
#       1: viewer_data: Path to qtl viewer .RData 
#       2: cc_variant:  Path to cc_variant.sqlite file
#       3: mgi_genes:   Path to mouse_genes_mgi.sqlite file
#       4: dataset:     Which dataset in the qtl viewer to use
#       5: type_expr:   Which expression dataset to use
#       6: type_peak:   Which lod.peaks data frame to use
#       7: perm_result: (Optional) Permutation Matrix as outputted by scan1perm. 'NA' or 'na' if not used.
#       8: alpha:       (Optional) but needed if given perm_result, see summary_scan1perm. Will subset peaks table to LOD scores above its alpha. 'NA' or 'na' if not used.
#       9: threshold:   (Optional) If given, this will subset peaks table based on LOD score. 'NA' or 'na' if not used.
#      10: pdf_file:    Name for the .pdf file to save
#
#
#   Output: 
#       1: Pdf file of all SNP Association plot for each QTL
#
#
#
#
#   Authors: Duy Pham
#   Date:    April 5, 2019
#   E-mails: duy.pham@jax.org
####################################################################################################################


### Options and Libraries
options(stringsAsFactors = FALSE)
library(qtl2)
library(dplyr)
library(qtl2utils)








### Load data
load('~/Desktop/Weinstock DOMA/weinstock_doma_16s_microbiome_qtl_viewer.RData')
cc_variant  <- '~/Desktop/cc_variants.sqlite'
mgi_genes   <- '~/Desktop/mouse_genes_mgi.sqlite'
dataset     <- 'dataset.doma.microbiome'
type_peak   <- 'additive'
type_expr   <- 'rankz'
perm_result <- '~/Desktop/Weinstock DOMA/weinstock_doma_16s_microbiome_additive_permutation.rds'
alpha       <- 0.05
threshold   <- 'NA'
pdf_file    <- '~/Desktop/weinstock_doma_snp_assocation_for_significant_peaks.pdf'








### Extract data
peaks <- get(dataset)$lod.peaks[[type_peak]]
expr  <- get(dataset)[[type_expr]]
covar <- get(dataset)$covar
query_variant <- create_variant_query_func(cc_variant)
query_genes   <- create_gene_query_func(mgi_genes)









### Keep LOD peaks above threshold if given perm_results
if(!(perm_result %in% c('NA','na'))){
   
   # Get permutation matrix and get LOD score at alpha level
   perm_result <- readRDS(perm_result)
   stopifnot(!(alpha %in% c('NA', 'na')))
   perm_summary <- summary_scan1perm(perm_result)
   
   
   
   # Find which peaks are above threshold and keep
   keep <- apply(X      = peaks, 
                 MARGIN = 1, 
                 FUN    = function(x) x['lod'] > perm_summary[,x['annot.id']])
  
   
   peaks <- peaks[keep, ]
}








### Keep peaks with LOD score above a given threshold
if(!(threshold %in% c('NA','na'))){
  
   peaks <- peaks[peaks$lod > threshold, ]
  
}










### SNP Association Plot
pdf(pdf_file, width = 11, height = 6)
for(i in 1:nrow(peaks)){
  
    snpasso <- scan1snps(genoprobs  = genoprobs,
                         map        = map,
                         pheno      = expr[,peaks$annot.id[i]],
                         kinship    = K,
                         addcovar   = dataset.doma.microbiome$covar,
                         query_func = query_variant,
                         chr        = peaks$qtl.chr[i],
                         start      = peaks$qtl.pos[i] - 5,
                         end        = peaks$qtl.pos[i] + 5)
  
    
    genes <- query_genes(chr = peaks$qtl.chr[i], start = peaks$qtl.pos[i] -5, end = peaks$qtl.pos[i] + 5)
    
    
    plot_snpasso2(scan1output = snpasso$lod,
                  snpinfo     = snpasso$snpinfo,
                  max_peak    = peaks$qtl.pos[i],
                  genes       = genes,
                  main        = paste0(peaks$annot.id[i],', Scan1 LOD = ', signif(peaks$lod[i])))
   
}

dev.off()






