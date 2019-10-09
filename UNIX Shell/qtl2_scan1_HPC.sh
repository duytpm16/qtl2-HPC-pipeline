####################################################################################################################
#
#   This script performs qtl2 additive or interactive scan using R script located:
#       https://github.com/duytpm16/qtl2-HPC-pipeline/blob/master/R%20scripts/qtl2_scan1.R
#
#   Notes*:
#         The QTL scan can be run in 'chunks' or all at once.
#             
#             Example for chunk size/number: 
#             Suppose there are 5433 phenotype columns. If chunk size is 1000, then there should be 6 different chunk_number,
#                  with the chunk_number value being 1,2,3,4,5, or 6 to get the columns:
#                       1-1000,1001-2000,2001-3000,3001-4000,4001-5000,5001-5433, respectively.
#
#         If you do not want to run in chunks, make chunk_number and chunk_size 'NA' or 'na'
#
#
#
#
#   Input (for shell script; see (https://github.com/duytpm16/qtl2-HPC-pipeline/blob/master/UNIX%20Shell/qtl2_scan1_HPC.sh):
#      1: viewer:        Path to the qtl viewer .RData 
#      2: dataset_expr:  Which dataset.* and data type to use. Ex. dataset.islet.proteins|data (if only one expression matrix) or dataset.islet.proteins|data|rz (if there are multiple expression matrices in 'data')
#      3: int_name:      (Optional) A string with interactive variables in samples dataframe separated by '|'.     Ex. 'sex' or 'sex|batch'. 'NA' or 'na' if not used
#      4: num_cores:     Number of cores to run
#      5: chunk_number:  (Optional) Numeric value of the chunk number. 'NA' or 'na' if not used
#      6: chunk_size:    (Optional) Numeric value of chunk size. Should be consistent. 'NA' or 'na' if not used
#
#   Output: 
#       1: Matrix containing LOD scores for each of the phenotype that was given to scan1 at each marker.
#
#
#
#
#   Authors: Duy Pham, Andrew Deighan, & Isabela Gyuricza
#   Date:    October 31, 2018
#   E-mails: duy.pham@jax.org, andrew.deighan@jax.org, & isabela.gyuricza@jax.org
####################################################################################################################


#PBS -q batch
#PBS -l nodes=1:ppn=1

module load R/3.5.1



job_name='esc_prot_sex'

for i in {1..79}
do
  echo "#PBS -l nodes=1:ppn=8
#PBS -q batch
#PBS -l walltime=72:00:00

module load R/3.5.1


viewer_data='munger_esc_proteomics_v2.RData'
dataset_expr='dataset.esc.proteins|data|rz'
num_cores='8'
int_name='sex'
chunk_number=${i}
chunk_size='100'


Rscript qtl2_scan1.R viewer_data=\$viewer_data dataset_expr=\$dataset_expr num_cores=\$num_cores int_name=\$int_name chunk_number=\$chunk_number chunk_size=\$chunk_size" >> "${job_name}_${i}.sh"
qsub "${job_name}_${i}.sh"
done
