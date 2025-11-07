#! /bin/bash
## pyscenic shell script

#load database
dir=/datf/mazhuo/database/scenic
tfs=$dir/allTFs_hg38.txt
feather=$dir/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
tbl=$dir/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl

celltypes=("Immature_AT1" "AT1" "AT2" "TRB_Secretory" "PreTB_Secretory" "Basal" "Goblet" "Ciliated")
for celltype in "${celltypes[@]}"
do
#input loom file
input_loom=/datf/mazhuo/jupyter_notebook/COPD/SCENIC/${celltype}_DEG_HCHS.loom
ls $tfs $feather $tbl

pyscenic grn \
        --num_workers 10 \
        --output /datf/mazhuo/jupyter_notebook/COPD/SCENIC/adj_${celltype}_HCHS.sample.tsv \
        --method grnboost2 \
        $input_loom $tfs

pyscenic ctx \
        /datf/mazhuo/jupyter_notebook/COPD/SCENIC/adj_${celltype}_HCHS.sample.tsv $feather \
        --annotations_fname $tbl \
        --expression_mtx_fname $input_loom \
        --mode "dask_multiprocessing" \
        --output /datf/mazhuo/jupyter_notebook/COPD/SCENIC/reg_${celltype}_HCHS.csv \
        --num_workers 10 \
        --mask_dropouts

pyscenic aucell \
        $input_loom \
        /datf/mazhuo/jupyter_notebook/COPD/SCENIC/reg_${celltype}_HCHS.csv \
        --output /datf/mazhuo/jupyter_notebook/COPD/SCENIC/out_SCENIC_${celltype}_DEG_HCHS.loom \
        --num_workers 10
done

