#PBS -q batch
#PBS -l nodes=1:ppn=1

module load R/3.5.1



job_name='doma_perm_age'

for i in {1..206}
do
  echo "#PBS -l nodes=1:ppn=4
#PBS -q batch
#PBS -l walltime=72:00:00

module load R/3.5.1


viewer_data='weinstock_doma_16s_microbiome_qtl_viewer.RData'
dataset='dataset.doma.microbiome'
num_cores='4'
type_data='rankz'
int_name='cohort.age'
perm_run='1000'
chunk_number=${i}
chunk_size='1'


Rscript qtl2_scan1perm.R viewer_data=\$viewer_data dataset=\$dataset num_cores=\$num_cores type_scan=\$type_scan type_data=\$type_data int_name=\$int_name perm_run=\$perm_run chunk_number=\$chunk_number chunk_size=\$chunk_size" >> "${job_name}_${i}.sh"
qsub "${job_name}_${i}.sh"
done
