import logging
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import numpy as np
import os
import scanpy as sc
from utils import check_matrix_files, output_adata


def qc(config):

    cache_dir = os.path.join(config['cache_dir'], 'cohort-level')
    input_path = os.path.join(cache_dir, 'merged.h5ad')

    if not check_matrix_files(input_path, 'qc'):
        logging.info('---------- Missing Input File: merged.h5ad ----------')        
        return
    
    figure_dir = os.path.join(config['figure_dir'], 'qc')
    os.makedirs(figure_dir, exist_ok=True)
    
    adata = sc.read_h5ad(input_path)
    logging.info('Merged adata shape: %s', adata.shape)

    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # ---------- Adjust QC metrics here ----------
    # ---------- TODO: Should be configured in config.json in the future ----------
    min_gene = config['qc_metrics']['min_gene_per_cell']
    max_doublet = config['qc_metrics']['max_doublet_score']
    max_mt = config['qc_metrics']['max_mt_percentage']

    adata.obs['Low_nFeature'] = adata.obs['n_genes_by_counts'] < min_gene
    adata.obs['Doublet'] = adata.obs['doublet_score'] > max_doublet
    adata.obs['High_MT'] = adata.obs['pct_counts_mt'] > max_mt
    adata.obs['Pass'] = (adata.obs['n_genes_by_counts'] >= min_gene) & (adata.obs['doublet_score'] <= max_doublet) & (adata.obs['pct_counts_mt'] <= max_mt)

    conditions = [
        adata.obs['Pass'] == True,
        adata.obs['Pass'] == False
    ]
    values = ['Pass', 'Failed']
    adata.obs['QC'] = np.select(conditions, values, default='Unknown')
    adata.obs['QC'] = adata.obs['QC'].astype('category')
    
    logging.info(adata.obs['QC'].value_counts())


    # Calculate pass ratio
    total_cells = len(adata.obs)
    pass_cells = sum(adata.obs['QC'] == 'Pass')
    pass_ratio = pass_cells / total_cells * 100

    # Get cell indices for each category
    doublet_cells = set(adata.obs[adata.obs['Doublet']].index)
    low_feature_cells = set(adata.obs[adata.obs['Low_nFeature']].index)
    high_mt_cells = set(adata.obs[adata.obs['High_MT']].index)

    # Create Venn diagram
    plt.figure(figsize=(10, 10))
    venn3([doublet_cells, low_feature_cells, high_mt_cells], 
        set_labels=('Doublet', 'Low_nFeature', 'High_MT'))
    plt.title(f'#cells failed any of the criteria\nPass Rate: {pass_ratio:.2f}%')
    plt.savefig(os.path.join(figure_dir, 'venn.png'))

    adata = adata[adata.obs['QC'] == 'Pass']

    output_adata(adata, os.path.join(cache_dir, 'qced.h5ad'))