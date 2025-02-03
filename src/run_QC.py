import logging
import numpy as np
import os
import scanpy as sc
from utils import check_matrix_files, output_adata


def qc(config):

    cache_dir = config['cache_dir']
    input_path = os.path.join(cache_dir, 'merged.h5ad')

    if not check_matrix_files(input_path, 'merge'):
        logging.info('---------- Missing Input File: merged.h5ad ----------')        
        return
    
    # output_dir = config['figure_dir']
    # os.makedirs(output_dir, exist_ok=True)
    # sc.settings.figdir = output_dir
    
    adata = sc.read_h5ad(input_path)
    logging.info('Merged adata shape: %s', adata.shape)

    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # ---------- Adjust QC metrics here ----------
    # ---------- TODO: Should be configured in config.json in the future ----------
    adata.obs['Low_nFeature'] = adata.obs['n_genes_by_counts'] <= 200
    adata.obs['Doublet'] = adata.obs['doublet_score'] >= 1.00
    adata.obs['High_MT'] = adata.obs['pct_counts_mt'] >= 1.00
    adata.obs['Pass'] = (adata.obs['n_genes_by_counts'] > 200) & (adata.obs['doublet_score'] < 1.00) & (adata.obs['pct_counts_mt'] < 1.00)
    # logging.info('Columns in adata.obs: %s', adata.obs.columns.tolist())

    conditions = [
        adata.obs['Low_nFeature'],
        adata.obs['Doublet'],
        adata.obs['High_MT'],
        adata.obs['Pass']
    ]
    values = ['Low_nFeature', 'Doublet', 'High_MT', 'Pass']
    adata.obs['QC'] = np.select(conditions, values, default='Unknown')
    adata.obs['QC'] = adata.obs['QC'].astype('category')
    
    logging.info(adata.obs['QC'].value_counts())

    adata = adata[adata.obs['QC'] == 'Pass']

    output_adata(adata, os.path.join(cache_dir, 'qced.h5ad'))