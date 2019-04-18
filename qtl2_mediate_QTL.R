### Options and Libraries
options(stringsAsFactors = FALSE)
library(dplyr)
library(tidyr)
library(intermediate)
library(intermediate2) # devtools::install_github('duytpm16/intermediate2')
library(qtl2)







### Command line arguments / variables to change. See .sh script too
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






### Load QTL viewer data. Changing some character parameters to numeric
load(viewer_data)
z_thres      <- as.numeric(z_thres)
pos_thres    <- as.numeric(pos_thres)
cores        <- as.numeric(cores)
targ_dataset_expr <- strsplit(targ_dataset_expr, '|',fixed = TRUE)[[1]]
med_dataset_expr  <- strsplit(med_dataset_expr, '|',fixed = TRUE)[[1]]
stopifnot(c('genoprobs', 'K', targ_dataset_expr[1], med_dataset_expr[1]) %in% ls())
targ_dataset <- targ_dataset_expr[1]
med_dataset  <- med_dataset_expr[1]





### Extract rankZ, annotation, covariate data, and target's lod.peaks table
targ_annot <- get(targ_dataset)[[targ_annot]] %>% dplyr::rename(id = targ_id, pos = start) %>% mutate(chr = as.character(chr))
med_annot  <- get(med_dataset)[[med_annot]] %>% dplyr::rename(id = med_id, pos = start) %>% mutate(chr = as.character(chr))



targ_expr <- get(targ_dataset)[[targ_dataset_expr[2]]]
med_expr  <- get(med_dataset)[[med_dataset_expr[2]]]

if(length(targ_dataset_expr) == 3){
   targ_expr <- targ_expr[[targ_dataset_expr[3]]][,targ_annot$id]
}

if(length(med_dataset_expr) == 3){
   med_expr <- med_expr[[med_dataset_expr[3]]][rownames(targ_expr),med_annot$id] 
}else{
   med_expr <- med_expr[rownames(targ_expr),med_annot$id] 
}


targ_covar <- get(targ_dataset)$covar.matrix[rownames(targ_expr),]
med_covar  <- get(med_dataset)$covar.matrix[rownames(targ_expr),]





lod.peaks  <- get(targ_dataset)$lod.peaks[[type_peak]]









### Making sure names match
stopifnot(colnames(targ_expr)  == targ_annot$id)
stopifnot(colnames(med_expr)   == med_annot$id)
stopifnot(rownames(targ_expr)  == rownames(med_expr))
stopifnot(rownames(med_expr)   == rownames(targ_covar))
stopifnot(rownames(targ_covar) == rownames(med_covar))












### Split lod.peaks table if runnning in parallel
targ_rng <- 1:nrow(lod.peaks)

if(!chunk_number %in% c('NA','na')){
  
   chunk_number <- as.numeric(chunk_number)
   chunk_size   <- as.numeric(chunk_size)
   
   max_col = nrow(lod.peaks)
   targ_rng = ((chunk_number - 1) * chunk_size + 1):(chunk_number * chunk_size)
   if(targ_rng[length(targ_rng)] > max_col) {
      targ_rng = targ_rng[1]:max_col
   }
}

lod.peaks <- lod.peaks[targ_rng,]












### Create vectors to store the data
#     Initialize empty vectors to store info
n <- nrow(lod.peaks)
results <- data.frame(target.id     = unname(lod.peaks[,targ_id]),
                      target.symbol = lod.peaks$gene.symbol,
                      target.chr    = lod.peaks$gene.chr,
                      target.start  = lod.peaks$gene.start,
                      target.end    = lod.peaks$gene.end,
                      qtl.chr       = lod.peaks$qtl.chr,
                      qtl.pos       = unname(lod.peaks[,qtl_pos_name]),
                      qtl.marker    = lod.peaks$marker.id,
                      mediator.id   = character(length = n),
                      mediator.symbol = character(length = n),
                      mediator.chr    = character(length = n),
                      mediator.start  = character(length = n),
                      mediator.end    = character(length = n),
                      pearson         = numeric(length = n),
                      best.model      = character(length = n),
                      best.model.p    = numeric(length = n),
                      best.triad      = character(length = n),
                      best.triad.p    = numeric(length = n),
                      causal.p        = numeric(length = n),
                      reactive.p      = numeric(length = n),
                      independent.p   = numeric(length = n),
                      undecided.p     = numeric(length = n),
                      target.lod      = lod.peaks$lod,
                      mediator.lod    = numeric(length = n),
                      mediation.lod   = numeric(length = n),
                      inverse.lod     = numeric(length = n),
                      mediation.z     = numeric(length = n),
                      inverse.z       = numeric(length = n))
















### For each QTL...
for(i in 1:nrow(results)){
    
  
    # Get distal QTL information
    target  <- results$target.id[i]
    qtl.chr <- results$qtl.chr[i]
    qtl.pos <- results$qtl.pos[i]
    marker  <- results$qtl.marker[i]
    
    
    
  
    
    
    
    
    
    # Do mediation scan on target at QTL against all mediators
    med <- mediation.scan(target     = targ_expr[, target, drop = FALSE],
                          mediator   = med_expr,
                          annotation = med_annot,
                          qtl.geno   = genoprobs[[qtl.chr]][rownames(targ_expr), , marker],
                          covar      = targ_covar,
                          verbose    = FALSE,
                          method     = med_method)
  
    
    
    
    
    
    
    
    
    
    # Filter mediation results that satisfy threshold requirement
    med <- med %>% 
               mutate(scaled_LOD = scale(LOD)) %>%
               filter(scaled_LOD < z_thres & chr == qtl.chr & abs(pos - qtl.pos) <= pos_thres)
    
  
    
    
  
    
    
    
    
    # If there are mediators that satisfy threshold...
    if(nrow(med) != 0){
        
       # Get mediator info
       results$mediator.id[i]     <- paste0(med$id,     collapse = ',')
       results$mediator.symbol[i] <- paste0(med$symbol, collapse = ',')
       results$mediator.chr[i]    <- paste0(med$chr,    collapse = ',')
       results$mediator.start[i]  <- paste0(med$pos,  collapse = ',')
       results$mediator.end[i]    <- paste0(med$end,    collapse = ',')
       results$mediation.lod[i]   <- paste0(signif(med$LOD),        collapse = ',')
       results$mediation.z[i]     <- paste0(signif(med$scaled_LOD), collapse = ',')
       
       
       
       
       
       
       # Compute correlation between target and mediator
       results$pearson[i]  <- paste0(signif(c(cor(targ_expr[,target], med_expr[,med$id], use = 'pairwise.complete.obs', method = 'pearson'))), collapse = ',')
      
  
    
    
    
       
       
    
       # Compute LOD score at QTL marker for mediators
       gp      <- genoprobs[,qtl.chr]
       gp[[1]] <- gp[[1]][rownames(targ_expr),,marker, drop = FALSE]  
       mediator.qtl <- scan1(genoprobs = gp, 
                             pheno     = med_expr[,med$id, drop = FALSE], 
                             kinship   = K[[qtl.chr]], 
                             addcovar  = med_covar, 
                             cores     = cores)
       results$mediator.lod[i] <- paste0(signif(c(mediator.qtl)), collapse = ',')
    

    
  
    
    
    
  
  
       # CausalMST and Inverse Mediation
       for(j in 1:nrow(med)){
         
           
           # Causal MST
           mst <- mediation_test(target     = targ_expr[,target,drop = FALSE],
                                 mediator   = med_expr[,med$id[j], drop = FALSE],
                                 driver     = genoprobs[[qtl.chr]][rownames(targ_expr), , marker],
                                 covar_tar  = targ_covar,
                                 covar_med  = med_covar)
          
           med$best.model[j]    <- as.character(mst$test$model[which.min(mst$test$pvalue)])
           med$best.model.p[j]  <- signif(as.numeric(mst$test$pvalue[which.min(mst$test$pvalue)]))
           med$best.triad[j]    <- as.character(mst$best$triad)
           med$best.triad.p[j]  <- signif(as.numeric(mst$best$pvalue))
           med$causal.p[j]      <- signif(as.numeric(mst$test[mst$test$model == 'causal',     'pvalue']))
           med$reactive.p[j]    <- signif(as.numeric(mst$test[mst$test$model == 'reactive',   'pvalue']))
           med$independent.p[j] <- signif(as.numeric(mst$test[mst$test$model == 'independent','pvalue']))
           med$undecided.p[j]   <- signif(as.numeric(mst$test[mst$test$model == 'undecided',  'pvalue']))
           
           
           
           
           
           # Inverse Mediation: M ~ Q + T + covar
           inv.med <- mediation.scan(target     = med_expr[, med$id[j], drop = FALSE],
                                     mediator   = targ_expr,
                                     annotation = targ_annot,
                                     qtl.geno   = genoprobs[[qtl.chr]][rownames(med_expr), , marker],
                                     covar      = med_covar,
                                     verbose    = FALSE,
                                     method     = med_method)
           inv.med <- inv.med %>% mutate(scaled_LOD = scale(LOD)) %>% filter(id == target)
           
           med$inverse.lod[j] <- signif(inv.med$LOD)
           med$inverse.z[j]   <- signif(inv.med$scaled_LOD)
          
        }

  
       
       
       
       
       
        # Store mediation results
        results$best.model[i]    <- paste0(med$best.model, collapse = ',')
        results$best.model.p[i]  <- paste0(med$best.model.p, collapse = ',')
        results$best.triad[i]    <- paste0(med$best.triad, collapse = ',')
        results$best.triad.p[i]  <- paste0(med$best.triad.p, collapse = ',')
        results$causal.p[i]      <- paste0(med$causal.p, collapse = ',')
        results$reactive.p[i]    <- paste0(med$reactive.p, collapse = ',')
        results$independent.p[i] <- paste0(med$independent.p, collapse = ',')
        results$undecided.p[i]   <- paste0(med$undecided.p, collapse = ',')
        results$inverse.lod[i]   <- paste0(med$inverse.lod, collapse = ',')
        results$inverse.z[i]     <- paste0(med$inverse.z, collapse = ',')
    }
    
    
    
    
    # Check 
    print(paste(i,'of', n))
}










### Make results long
results <- results %>%
                   filter(mediator.id != '') %>%
                   separate_rows(mediator.id, mediator.symbol, mediator.chr, mediator.start, mediator.end, 
                                 pearson, best.model, best.model.p, best.triad, best.triad.p,
                                 causal.p, reactive.p, independent.p, undecided.p,
                                 mediator.lod, mediation.lod, inverse.lod, mediation.z,
                                 inverse.z, sep = ',') %>%
                   mutate(mediator.start = as.numeric(mediator.start),
                          mediator.end   = as.numeric(mediator.end),
                          pearson        = as.numeric(pearson),
                          best.model.p   = as.numeric(best.model.p),
                          best.triad.p   = as.numeric(best.triad.p),
                          causal.p       = as.numeric(causal.p),
                          reactive.p     = as.numeric(reactive.p),
                          independent.p  = as.numeric(independent.p),
                          undecided.p    = as.numeric(undecided.p),
                          mediator.lod   = as.numeric(mediator.lod),
                          mediation.lod  = as.numeric(mediation.lod),
                          inverse.lod    = as.numeric(inverse.lod),
                          mediation.z    = as.numeric(mediation.z),
                          inverse.z      = as.numeric(inverse.z))







### Add LOD drop proportion
results$mediation.dp <- (results$target.lod - results$mediation.lod) / results$target.lod
results$inverse.dp   <- (results$mediator.lod - results$inverse.lod) / results$mediator.lod













### Save data as .rds file
if(chunk_size != 'NA'){
   saveRDS(results, file = paste0(filename,'_chunk_',chunk_number,'.rds'))
}else{
   saveRDS(results, file = paste0(filename,'.rds'))
}
