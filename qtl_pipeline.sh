#PBS -l nodes=1:ppn=8
#PBS -q batch
#PBS -l walltime=72:00:00


module load R/3.5.1


### Uncomment line 10 to run scan1 not in parallel, then comment line 11. Otherwise qsub run_scan1_parallel.sh first,
#       then qsub this .sh file
#Rscript scan1_HPC.R attie_islet_284_qtl_viewer_v2.RData dataset.islet.mrna rankz 8
Rscript gather_scan1_chunks.R islet_mrna_additive_scan 1 22 test.rds
Rscript find_peaks.R attie_islet_284_qtl_viewer_v2.RData attie_islet_mrna_rZ_qtl_lod.rds dataset.islet.mrna 6 8 additive rankz
