#PBS -l nodes=1:ppn=8
#PBS -q batch
#PBS -l walltime=72:00:00
module load R/3.5.1


pattern='esc_proteins_additive_scan'
chunk_start='1'
chunk_end='197'
func='cbind'
out_file='munger_esc_proteins_additive_qtl_lod'



Rscript gather_scan1_chunks.R "pattern=$pattern" "chunk_start=$chunk_start" "chunk_end=$chunk_end" "func=$func" "out_file=$out_file"
