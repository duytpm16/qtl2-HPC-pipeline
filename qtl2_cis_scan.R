options(stringsAsFactors = FALSE)
library(qtl2)



### Variables to change
#    * - https://github.com/churchill-lab/qtlapi/blob/master/docs/QTLAPIDataStructures.md
#    Assuming genoprobs, map, and markers are in the environment

#    1.) annots   : annotation dataframe as defined in *
#    2.) covar    : matrix of covariates as created by model.matrix.  samples x covariates
#    3.) expr     : matrix of expressions.                            samples x expressions
#    4.) id_col   : column name in annots dataframe that correspond to the column name of expr
#    5.) chr_col  : name of column in annots that correspond to the chromosome location of each expr
#    6.) pos_col  : name of column in annots that correspond to the starting Mb location of each expr.
#    7.) cis_thres: maximum distance from pos_col to determine cis.  cis = (chr_col == marker$chr & abs(pos_col - marker$pos) <= cis_thres)
#    8.) outfile  : name to save results as .rds
annots    <- dataset.islet.mrna$annots
covar     <- dataset.islet.mrna$covar
expr      <- dataset.islet.mrna$rankz
id_col    <- 'gene_id'
chr_col   <- 'chr'
pos_col   <- 'start'
cis_thres <- 4
outfile   <- 'attie_islet_mrna_cis_scan_results.rds'




### Extracting required values from annots 
ids    <- annots[,id_col]       
id_chr <- annots[,chr_col]
id_pos <- annots[,pos_col]
     




### Creating dataframe to save results
n <- nrow(annots) 
results <- data.frame(id_col      = ids,
                      chr_col     = id_chr, 
                      pos_col     = id_pos,
                      best_marker = character(length = n),
                      lod         = numeric(length = n))
                      
results <- cbind(results, 
                 matrix(0, nrow = nrow(annots), ncol = 8, dimnames = list(1:nrow(annots), LETTERS[1:8])))







### cis Scan begins
for(i in 1:nrow(annots)){
       
    # Find markers within cis_threshold to the start position of an expression 
    nearest_markers <- subset(markers, chr == id_chr[i] & abs(pos-id_pos[i] <= cis_thres))$marker
        
        
    
    # Scan1 on all markers closest to the start position of an expression
    gp        <- genoprobs[,id_chr[i]]
    gp[[1]]   <- gp[[1]][,,nearest_markers, drop = FALSE]
    cis_scan1 <- scan1(genoprobs = gp, 
                       pheno     = expr[,ids[i], drop = FALSE],
                       kinship   = K[[id_chr[i]]],
                       addcovar  = covar)
        
    
    
    # Get cis marker with highest LOD score
    best_cis_marker <- max_scan1(scan1_output = cis_scan1, 
                                 map          = map)
    
    
    
    # Get allele effect at best cis marker
    gp[[1]] <- gp[[1]][,,rownames(best_cis_marker),drop = FALSE]
    blup    <- scan1blup(genoprobs = gp,
                         pheno     = expr[,ids[i], drop = FALSE],
                         kinship   = K[[id_chr[i]]],
                         addcovar  = covar)
    
    
    
    
    
    # Save LOD, marker name, and allele effects at the best cis marker
    results$lod[i]           <- best_cis_marker[1, ids[i]]
    results$best_marker[i]   <- rownames(best_cis_marker) 
    results[i, LETTERS[1:8]] <- blup[1, LETTERS[1:8]]
  
    
    print(i)
  
}




colnames(results)[1:3] <- c(id_col, chr_col, pos_col)
saveRDS(results, file = outfile)
