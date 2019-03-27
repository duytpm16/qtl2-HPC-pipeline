#PBS -l nodes=1:ppn=8
#PBS -l walltime=24:00:00


module load R/3.5.1


job_name='attie_distal_pQTL_mediation'

for i in {1..68}
do
  echo "#PBS -l nodes=1:ppn=8
#PBS -q batch
#PBS -l walltime=24:00:00

module load R/3.5.1



viewer_data='attie_islet_284_qtl_viewer_v2.RData'
targ_dataset='dataset.islet.proteins.284'
med_dataset='dataset.islet.mrna.284'
targ_id='protein_id'
med_id='gene_id'
expr_type='rankz'
med_method='double-lod-diff'
z_thres='-4'
pos_thres='10'
cores='8'
filename='attie_islet_proteins_distal_pQTL_mediation'
chunk_number=${i}
chunk_size='100'






Rscript mediation_distal_QTL.R viewer_data=\$viewer_data targ_dataset=\$targ_dataset med_dataset=\$med_dataset targ_id=\$targ_id med_id=\$med_id expr_type=\$expr_type med_method=\$med_method z_thres=\$z_thres pos_thres=\$pos_thres cores=\$cores filename=\$filename chunk_number=\$chunk_number chunk_size=\$chunk_size" >> "${job_name}_${i}.sh"
qsub "${job_name}_${i}.sh"
done




