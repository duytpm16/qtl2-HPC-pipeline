#PBS -q batch
#PBS -l nodes=1:ppn=1

module load R/3.5.1



job_name='doma_perm'

for i in {1..21}
do
  echo "#PBS -l nodes=1:ppn=8
#PBS -q batch
#PBS -l walltime=72:00:00

module load R/3.5.1

viewer_data='weinstock_doma_16s_microbiome_qtl_viewer.RData'
dataset_expr='dataset.doma.microbiome|data|rz'
num_cores='8'
int_name='NA'
perm_run='1000'
chunk_number=${i}
chunk_size='10'

Rscript qtl2_scan1perm.R viewer_data=\$viewer_data dataset_expr=\$dataset_expr num_cores=\$num_cores int_name=\$int_name perm_run=\$perm_run chunk_number=\$chunk_number chunk_size=\$chunk_size" >> "${job_name}_${i}.sh"
qsub "${job_name}_${i}.sh"
done
