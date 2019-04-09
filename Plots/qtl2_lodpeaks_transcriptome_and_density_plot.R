####################################################################################################################
#
#   This script is used to generate a 3-layer transcriptome, cis QTL density, and distal QTL density plots.
#
#
#
#   Input:
#       1: viewer_data: Path to qtl viewer .RData 
#       2: lod.peaks:   Path to cc_variant.sqlite file
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
library(tidyverse)
library(grid)
library(gridExtra)
library(plyr)
library(data.table)





### Variables to change
viewer_data <- "munger_esc_proteomics_qtl_viewer.RData"
dataset     <- 'dataset.esc.proteins'
type_peak   <- 'additive'
slide  <- 1
window <- 4
lod.thres <- 7
density.thres <- 50
cis_color <- "red"
dis_color <- 'blue'






### Get and filter lod.peaks table
lod.peaks    <- get(dataset)[['lod.peaks']][[type_peak]]
lod.peaks    <- lod.peaks[lod.peaks$lod > lod.thres,]
lod.peaks    <- lod.peaks[complete.cases(lod.peaks),]
lod.peaks    <- lod.peaks[lod.peaks$gene.chr %in% c(1:19,'X'),]
local.peaks  <- lod.peaks[lod.peaks$cis,]
distal.peaks <- lod.peaks[!lod.peaks$cis,]
colnames(lod.peaks) <- gsub('.','_',colnames(lod.peaks), fixed = TRUE)












### Format chromosomes for transcriptome plot
all.chr = lod.peaks %>%
                    select(qtl_chr, gene_chr) %>%
                    gather(k, v) %>%
                    select(v) %>%
                    distinct() %>%
                    arrange(v)

all.chr <- all.chr$v[!is.na(all.chr$v)]

if(length(grep("M", all.chr)) > 0){
   wh <- grep("M", all.chr)
   all.chr <- all.chr[c(1:(wh-1), (wh+1):length(all.chr), wh)]
}




### Format lod.peaks for transcriptome plot
data = lod.peaks %>% 
                 mutate(cis      = (gene_chr == qtl_chr) & (abs(gene_start - qtl_pos) <= 4),
                        qtl_chr  = factor(qtl_chr, levels   = all.chr[order(as.numeric(all.chr))]),
                        gene_chr = factor(gene_chr,levels   = rev(all.chr[order(as.numeric(all.chr))])),
                        gene_pos = (gene_end + gene_start) * 0.5)
cis.colors        = c(dis_color, cis_color)
names(cis.colors) = c("FALSE", "TRUE")















## Cis LOD counts
lod_df_cis <- list()
for(i in unique(markers$chr)){
  
    # Finding floor of minimum marker position and ceiling of maximum marker position
    min <- round_any(min(map[[i]]), 1, f = floor)
    max <- round_any(max(map[[i]]), 4, f = ceiling)
    
    # Creating x-axis scale. min to max with slide (or 'by' in seq function)
    x_axis_scale <- seq(min, max, slide)
    chr <- rep(i, length(x_axis_scale))
    
    # Getting LOD peaks from chromosome i
    sub <- subset(local.peaks, qtl.chr == i)
    
    
    # Creating dataframe of counts of lod peaks above threshold 
    count <- vector()
    pos <- ((x_axis_scale+window)+x_axis_scale)/2
    for(j in 1:length(pos)){
        count[j] <- sum((abs(sub$qtl.pos - pos[j])  <= window/2))
    }
    
    lod_df_cis[[i]] <- data.frame(chr = chr, pos = pos, count = count)
}

lod_df_cis     <- rbindlist(lod_df_cis)
lod_df_cis$chr <- factor(lod_df_cis$chr, levels = c(1:19,'X'))











### Distal LOD counts
lod_df_dis <- list()
for(i in unique(markers$chr)){
  
    # Finding floor of minimum marker position and ceiling of maximum marker position
    min <- round_any(min(map[[i]]), 1, f = floor)
    max <- round_any(max(map[[i]]), 4, f = ceiling)
    
    # Creating x-axis scale. min to max with slide (or 'by' in seq function)
    x_axis_scale <- seq(min, max, slide)
    chr <- rep(i, length(x_axis_scale))
    
    # Getting LOD peaks from chromosome i
    sub <- subset(distal.peaks, qtl.chr == i)
    
    
    # Creating dataframe of counts of lod peaks above threshold 
    count <- vector()
    pos <- ((x_axis_scale+window)+x_axis_scale)/2
    for(j in 1:length(pos)){
        count[j] <- sum((abs(sub$qtl.pos - pos[j])  <= window/2))
    }
    
    lod_df_dis[[i]] <- data.frame(chr = chr, pos = pos, count = count)
  }

lod_df_dis <- rbindlist(lod_df_dis)
lod_df_dis$chr <- factor(lod_df_dis$chr, levels = c(1:19,'X'))


















### Plot cis density plot
c <- ggplot(lod_df_cis, aes(x = pos, y = count)) +
        geom_line(col = cis_color) +
        geom_hline(yintercept = density.thres, linetype = 2, color = 'grey70') +
        ylab(paste('No. localQTL/',window,'Mbp')) +
        xlab('Chromosome') +
        theme(panel.spacing.x = unit(0, "lines"),
              plot.margin = unit(c(0,0,1.3,0),"cm"),
              panel.background = element_blank(),
              panel.border = element_rect(fill = 0, color = "grey70"),
              axis.text.x = element_blank(),
              strip.text = element_text(size= 23),
              axis.ticks.x = element_blank(),
              strip.background = element_blank(),
              strip.text.x = element_text(size = 12)) + 
        scale_y_continuous(limits = c(0,80), breaks = seq(0,80,20)) +
        facet_grid(.~chr, scales = "free",switch = 'x')








### Plot distal density plot
d <- ggplot(lod_df_dis, aes(x = pos, y = count)) +
        geom_line(col = dis_color) +
        geom_hline(yintercept = density.thres, linetype = 2, color = 'grey70') +
        ylab(paste('No. distalQTL/',window,'Mbp')) +
        theme(panel.spacing.x = unit(0, "lines"),
              plot.margin = unit(c(2,0,0,0),"cm"),
              panel.background = element_blank(),
              panel.border = element_rect(fill = 0, color = "grey70"),
              axis.text.x = element_blank(),
              strip.text = element_text(size= 23),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              strip.background = element_blank(),
              strip.text.y = element_blank(),
              strip.text.x = element_blank()) +
        scale_y_continuous(limits = c(0,80), breaks = seq(0,80,20)) +
        facet_grid(.~chr, scales = "free")










### Plot transcriptome
g <- ggplot(data, aes(x = qtl_pos, y = gene_pos, size = log(lod)), alpha = 0.4) +
          geom_point(aes(color = cis), alpha = 0.4) + 
          scale_color_manual(values = cis.colors, labels = c('Distal-pQTL','Local-pQTL')) +
          facet_grid(gene_chr ~ qtl_chr, scales = "free",switch = 'y', shrink = TRUE) +
          theme(panel.background = element_blank(),
                panel.border = element_rect(fill = 0, color = "grey70"),
                panel.grid.minor = element_blank(),
                panel.spacing = unit(0.0, "line"),
                axis.ticks = element_blank(),
                strip.background = element_blank(),
                strip.text.x = element_blank(),
                strip.text.y = element_text(angle = 180),
                axis.text = element_blank(),
                axis.line.y.left = element_blank(),
                legend.position = 'none',
                legend.title = element_blank(),
                legend.background = element_blank(),
                legend.text = element_blank(),
                legend.key = element_blank()) +
          xlab('QTL position (Chr)') +
          ylab('Gene position (Chr)')















### Put plots in 3 panels vertically. DONT LOAD COWPLOT LIBRARY!!
cowplot::plot_grid(g, d,c, axis = 'lr', ncol = 1, align = 'v', rel_heights = c(1,.3,.3), labels = c('A','B','C'), vjust = c(1.5,6.0,.85))

grid.text(label = 'Distal-QTL', x=unit(.15,'npc'), y = unit(.97,'npc'),gp=gpar(fontsize=20, fontface = 'bold', col=dis_color))
grid.text(label = 'Local-QTL', x=unit(.15,'npc'), y = unit(.95,'npc'), gp=gpar(fontsize=20, fontface = 'bold', col=cis_color))
