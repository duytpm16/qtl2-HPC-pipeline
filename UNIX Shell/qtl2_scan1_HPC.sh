#PBS -l nodes=1:ppn=8
#PBS -q batch

module load R/3.5.1



job_name='run_esc_prot'

for i in {1..197}
do
  echo "#PBS -l nodes=1:ppn=8
#PBS -q batch

module load R/3.5.1


viewer_data='munger_esc_viewer_174_v3.RData'
dataset='dataset.esc.proteins'
num_cores='8'
type_data='rankz'
int_name='NA'
chunk_number=${i}
chunk_size=40


Rscript qtl2_scan1.R viewer_data=\$viewer_data dataset=\$dataset num_cores=\$num_cores type_scan=\$type_scan type_data=\$type_data int_name=\$int_name chunk_number=\$chunk_number chunk_size=\$chunk_size" >> "${job_name}_${i}.sh"
qsub "${job_name}_${i}.sh"
done
