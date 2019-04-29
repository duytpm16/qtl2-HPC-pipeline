####################################################################################################################
#
#   This script runs the find_peaks function in qtl2.
#
#
#
#   Input:
#       1: viewer_data:  Path to qtl viewer .RData 
#       2: scan1:        Path to scan1 matrix outputted from scan1 function as .rds (LOD matrix)
#       3: dataset_expr: Which dataset.* and data type to use. Ex. dataset.islet.proteins|data (if only one expression matrix) or dataset.islet.proteins|data|rz (if there are multiple expression matrices in 'data')
#       4: thr:          See thre parameter in find_peaks function.
#       5: num_cores:    Number of cores to run
#       6: type_scan:    Type of scan. Ex. 'additive', 'sex_int', 'age_int'...
#       7: drop:         (Optional) See 'drop' parameter of find_peaks function. Leave as 'NA' or 'na' if not used
#	8: int_mat:      (Optional) LOD matrix from interaction scan to get effects. Leave as 'NA' or 'na' if not used
#
#
#   Output: 
#       1: Output from find peaks function saved in the dataset given
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
print(args)
if(length(args)==0){
    print("No arguments supplied.")
}else{
    for(i in 1:length(args)){
        a <- strsplit(args[i],split = '=', fixed = TRUE)[[1]]
	assign(a[1],a[2])
    }
}








### Get required data
load(viewer_data)


dataset_expr <- strsplit(dataset_expr, '|', fixed = TRUE)[[1]]

stopifnot(c('map', dataset_expr[1]) %in% ls())

ds    <- get(dataset_expr[1])
expr  <- ds[[dataset_expr[2]]]

if(length(dataset_expr) == 3){
   expr <- expr[[dataset_expr[3]]]
}



scan1_mat <- readRDS(scan1_mat)
thr       <- as.numeric(thr)
num_cores <- as.numeric(num_cores)
cis_threshold <- as.numeric(cis_threshold)










### Get confidence intervals for find_peaks
if(!drop %in% c('NA','na')){
    drop <- as.numeric(drop)
}else{
    drop <- NULL
}












### Get interaction matrix if type_scan is not additive
if(!int_mat %in% c('NA','na')){
   int_mat <- readRDS(int_mat)
   
   stopifnot(dim(int_mat) == dim(scan1_mat))
   stopifnot(colnames(int_mat) == colnames(scan1_mat))
   stopifnot(rownames(int_mat) == rownames(scan1_mat))   

   scan1_mat <- int_mat - scan1_mat
   print(sum(scan1_mat <= 0))
   #stopifnot(sum(scan1_mat <= 0, na.rm=TRUE) == 0) 
}









### Run fin_peaks
peaks <- find_peaks(scan1_mat, 
		    map = map, 
		    threshold = thr, 
                    drop = drop, 
                    cores = num_cores)












### Formatting for QTL Viewer
if(!is.null(drop)){

   peaks <- peaks %>%
                  select(-lodindex) %>%                         # Remove lodindex column
                  left_join(markers, by = c('chr', 'pos')) %>%  # Get marker ids
                  rename(annot.id  = lodcolumn,
                         qtl.chr   = chr,
                         qtl.pos   = pos,
                         marker.id = marker,
			 ci.lo     = ci_lo,
			 ci.hi     = ci_hi) %>%
               	  select(annot.id, marker.id, lod, qtl.chr, qtl.pos, ci.lo, ci.hi)       # Reorder data frame
}else{

   peaks <- peaks %>%
                  select(-lodindex) %>%                         # Remove lodindex column
                  left_join(markers, by = c('chr', 'pos')) %>%  # Get marker ids
                  rename(annot.id  = lodcolumn,
                         qtl.chr   = chr,
                         qtl.pos   = pos,
                         marker.id = marker) %>%
                  select(annot.id, marker.id, lod, qtl.chr, qtl.pos)       # Reorder data frame


}











### Check if all marker id are present in markers data frame
#  If not, replace with closest marker
if(!all(peaks$marker.id %in% markers$marker)){
	
   index <- which(!peaks$marker.id %in% markers$marker)
   
   for(i in index){

       sub_markers <- subset(markers, chr == peaks$qtl.chr[i])
       sub_markers <- sub_markers[which.min(abs(sub_markers$pos - peaks$qtl.pos[i])),]
	
       peaks$qtl.chr[i] <- sub_markers$chr
       peaks$qtl.pos[i] <- sub_markers$pos
       peaks$marker.id[i] <- sub_markers$marker
	   
   }

}












### BLUP scan at each QTL if additive
if(type_scan == 'additive'){
   stopifnot(all(peaks$marker.id %in% markers$marker))
   stopifnot(c('K','genoprobs') %in% ls())
  
   peaks = cbind(peaks, matrix(0, nrow = nrow(peaks), ncol = 8,
                               dimnames = list(NULL, LETTERS[1:8])))
   covar = ds$covar.matrix
   

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










### Determine if QTLs are cis or not
if((ds$datatype %in% c('mRNA', 'protein'))){

   id = if(ds$datatype == 'mRNA') 'gene_id' else 'protein_id'
   annots = if(ds$datatype == 'mRNA') 'annot.mrna' else 'annot.protein'
   annots <- ds[[annots]]
   peaks  <- merge(peaks, annots[,c(id,'chr','start','end','symbol')], by.x='annot.id',by.y = id, all.x = TRUE)

  
   peaks <- peaks %>% 
		  rename(gene.start=start,
			 gene.end=end,
		         gene.chr=chr,
			 gene.symbol=symbol) %>% 
		  mutate(gene.middle = (gene.start + gene.end)/2, cis = (gene.chr == qtl.chr) & abs(qtl.pos-gene.middle) <= cis_threshold)

   print(head(peaks))

   if(!is.null(drop)){
      
      if(is.matrix(int_mat)){

      	 peaks <- peaks %>% 
	         	select(annot.id, marker.id, lod, qtl.chr, qtl.pos, gene.chr, gene.start, gene.end, gene.symbol, cis, ci.lo, ci.hi)
      }else{ 
         
        peaks <- peaks %>%
                        select(annot.id, marker.id, lod, qtl.chr, qtl.pos, gene.chr, gene.start, gene.end, gene.symbol, cis, ci.lo, ci.hi, A, B, C, D, E, F, G, H) 
         
      }
		  
   }else{

      if(is.matrix(int_mat)){

         peaks <- peaks %>%
                        select(annot.id, marker.id, lod, qtl.chr, qtl.pos, gene.chr, gene.start, gene.end, gene.symbol, cis)
      }else{

        peaks <- peaks %>%
                        select(annot.id, marker.id, lod, qtl.chr, qtl.pos, gene.chr, gene.start, gene.end, gene.symbol, cis, A, B, C, D, E, F, G, H)

      }



   }
}







annot.id <- switch(ds$datatype,
		   'protein'   = 'protein.id',
		   'mrna'      = 'gene.id',
		   'phenotype' = 'data.name')
colnames(peaks)[1] <- as_tibble(annot.id)








### Save peaks to data
ds$lod.peaks[[type_scan]] <- peaks
assign(dataset_expr[1], ds)










### Save to .RData file
save(list = ls()[ls() %in% c(grep('dataset[.]', ls(), value = TRUE),
                                  'K',
                                  'genoprobs',
                                  'map',
                                  'markers',
                                  'ensemble.version')], file = viewer_data)
