#PBS -l nodes=1:ppn=8
#PBS -q batch

module load R/3.5.1


viewer_data='prelim_proteomics_correctedIDs_v3.1.RData'
scan1_mat='munger_esc_proteins_additive_qtl_lod.rds'
dataset='dataset.esc.proteins'
thr='6'
num_cores='8'
type_scan='additive'
type_data='rankz'
drop='NA'
int_mat='NA'



Rscript qtl2_findpeaks.R "viewer_data=$viewer_data" "scan1_mat=$scan1_mat" "dataset=$dataset" "thr=$thr" "num_cores=$num_cores" "type_scan=$type_scan" "type_data=$type_data drop=$drop int_mat=$int_mat"

