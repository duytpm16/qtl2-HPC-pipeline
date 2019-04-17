#PBS -l nodes=1:ppn=8
#PBS -q batch
#PBS -l walltime=72:00:00
module load R/3.5.1


pattern='doma_microbiome_additive_scan_chunk'
chunk_start='1'
chunk_end='21'
func='cbind'
out_file='weinstock_doma_16s_microbiome_additive_qtl_lod'



Rscript qtl2_gather_chunks.R "pattern=$pattern" "chunk_start=$chunk_start" "chunk_end=$chunk_end" "func=$func" "out_file=$out_file"

