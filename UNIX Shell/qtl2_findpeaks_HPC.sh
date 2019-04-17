#PBS -l nodes=1:ppn=32
#PBS -q batch
#PBS -l walltime=72:00:00
module load R/3.5.1


viewer_data='weinstock_doma_16s_microbiome_qtl_viewer.RData'
scan1_mat='weinstock_doma_16s_microbiome_additive_qtl_lod.rds'
dataset_expr='dataset.doma.microbiome|data|rz'
thr='6'
num_cores='32'
type_scan='sex'
drop='1.5'
int_mat='weinstock_doma_16s_microbiome_sex_int_qtl_lod.rds'



Rscript qtl2_findpeaks.R "viewer_data=$viewer_data" "scan1_mat=$scan1_mat" "dataset_expr=$dataset_expr" "thr=$thr" "num_cores=$num_cores" "type_scan=$type_scan" "drop=$drop" "int_mat=$int_mat"

