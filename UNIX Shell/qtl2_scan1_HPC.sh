#PBS -q batch
#PBS -l nodes=1:ppn=1

module load R/3.5.1



job_name='esc_prot_sex'

for i in {1..79}
do
  echo "#PBS -l nodes=1:ppn=8
#PBS -q batch

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
