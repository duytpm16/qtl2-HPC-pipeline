for i in {1..22}
do
  echo "PBS -q batch
#PBS -l nodes=1:ppn=8

module load R/3.5.1

Rscript scan1_HPC.R attie_islet_284_qtl_viewer_v2.RData dataset.islet.mrna rankz 8  $i 1000">> run_mrna_$i.sh
qsub run_mrna_$i.sh
done
