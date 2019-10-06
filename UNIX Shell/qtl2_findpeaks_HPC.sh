####################################################################################################################
#
#   This shell script runs the find_peaks function in qtl2 on HPC.
#
#
#
#   Input:
#       1: viewer_data:  Path to qtl viewer .RData 
#       2: scan1:        Path to additive scan1 matrix. Outputted from scan1 function and saved as .rds file.
#       3: dataset_expr: Which dataset.* and data type to use. Ex. dataset.islet.proteins|data (if only one expression matrix) or dataset.islet.proteins|data|rz (if there are multiple expression matrices in 'data')
#       4: thr:          See thres parameter in find_peaks function.
#       5: num_cores:    Number of cores to run
#       6: type_scan:    Name to save peaks table in lod.peaks list. Type of scan. Ex. 'additive', 'sex_int', 'age_int'...
#       7: prob:         (Optional) See 'prob' parameter of find_peaks function. Leave as 'NA' or 'na' if not used
#	      8: int_mat:      (Optional) LOD matrix from interaction scan to get effects. Leave as 'NA' or 'na' if not used
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



#PBS -l nodes=1:ppn=32
#PBS -q batch
#PBS -l walltime=72:00:00
module load R/3.5.1


viewer_data='weinstock_doma_16s_microbiome_qtl_viewer.RData'
scan1_mat='weinstock_doma_16s_microbiome_additive_qtl_lod.rds'
dataset_expr='dataset.doma.microbiome|data|rz'
thr='6'
num_cores='32'
type_scan='sex_int'
drop='1.5'
cis_threshold='4'
int_mat='weinstock_doma_16s_microbiome_sex_int_qtl_lod.rds'



Rscript qtl2_findpeaks.R "viewer_data=$viewer_data" "scan1_mat=$scan1_mat" "dataset_expr=$dataset_expr" "thr=$thr" "num_cores=$num_cores" "type_scan=$type_scan" "drop=$drop" "cis_threshold=$cis_threshold" "int_mat=$int_mat"

