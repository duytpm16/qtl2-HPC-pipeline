#PBS -q batch
#PBS -l nodes=1:ppn=1

module load R/3.5.1



job_name='esc_prot_run'

for i in {1..80}
do
  echo "#PBS -l nodes=1:ppn=8
#PBS -q batch

module load R/3.5.1


viewer_data='munger_esc_proteomics_qtl_viewer_v1.RData'
dataset='dataset.esc.proteins'
num_cores='8'
type_data='rankz'
int_name='NA'
chunk_number=${i}
chunk_size='100'


Rscript qtl2_scan1.R viewer_data=\$viewer_data dataset=\$dataset num_cores=\$num_cores type_scan=\$type_scan type_data=\$type_data int_name=\$int_name chunk_number=\$chunk_number chunk_size=\$chunk_size" >> "${job_name}_${i}.sh"
qsub "${job_name}_${i}.sh"
done
