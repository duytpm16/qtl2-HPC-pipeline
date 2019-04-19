# viewer_data:  path to QTL formatted .RData
# targ_dataset: dataset.* of the target
# med_dataset:  dataset.* of the mediator
# targ_id:      target identification in annots data frame: Ex. 'protein_id' or 'gene_id'
# med_id:       mediator identification in annots data frame: Ex. 'protein_id' or 'gene_id'
# expr_type:    Which expression matrix to use in the dataset.*
# type_peak:    Which LOD summary table to use in lod.peaks list in dataset.*
# med_method:   See 'method' parameter in intermediate package by Petr Simecek
# z_thres:      Z-score cut-off. Remove mediators that do not drop LOD score below z_thres. Should be negative value
# pos_thes:     Keep mediators that are within some Mbp of the QTL
# cores:        Number of cores
# filename:     Name to save results
# chunk_size:   Chunk size. Should be consistent.
# chunk_number: Chunk number. Changes with 'i' below




#PBS -l nodes=1:ppn=8
#PBS -l walltime=24:00:00


module load R/3.5.1


job_name='attie_islet_pQTL_mediation'

for i in {1..79}
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
type_peak='additive'
med_method='double-lod-diff'
z_thres='-4'
pos_thres='10'
cores='8'
filename='attie_islet_pQTL_mediation'
chunk_number=${i}
chunk_size='100'
Rscript qtl2_mediate_QTL.R viewer_data=\$viewer_data targ_dataset=\$targ_dataset med_dataset=\$med_dataset targ_id=\$targ_id med_id=\$med_id expr_type=\$expr_type type_peak=\$type_peak med_method=\$med_method z_thres=\$z_thres pos_thres=\$pos_thres cores=\$cores filename=\$filename chunk_number=\$chunk_number chunk_size=\$chunk_size" >> "${job_name}_${i}.sh"
qsub "${job_name}_${i}.sh"
done



