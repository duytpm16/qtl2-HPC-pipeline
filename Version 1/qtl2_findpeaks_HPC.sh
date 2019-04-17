#PBS -l nodes=1:ppn=32
#PBS -q batch
#PBS -l walltime=72:00:00
module load R/3.5.1


viewer_data='munger_esc_proteomics_qtl_viewer_v1.RData'
scan1_mat='munger_esc_proteins_additive_qtl_lod_v1.rds'
dataset='dataset.esc.proteins'
thr='6'
num_cores='32'
type_scan='additive'
type_data='rankz'
drop='1.5'
int_mat='NA'



Rscript qtl2_findpeaks.R "viewer_data=$viewer_data" "scan1_mat=$scan1_mat" "dataset=$dataset" "thr=$thr" "num_cores=$num_cores" "type_scan=$type_scan" "type_data=$type_data" "drop=$drop" "int_mat=$int_mat"


