#PBS -l nodes=1:ppn=8
#PBS -q batch

module load R/3.5.1


viewer_data='prelim_proteomics_correctedIDs_v3.1.RData'
dataset='dataset.esc.proteins'
num_cores='8'
type_data='rankz'
int_name='sex'
chunk_number='NA'
chunk_size='NA'



Rscript qtl2_interaction_scan.R "viewer_data=$viewer_data" "dataset=$dataset" "num_cores=$num_cores" "type_data=$type_data" "int_name=$int_name" "chunk_number=$chunk_number" "chunk_size=$chunk_size"
