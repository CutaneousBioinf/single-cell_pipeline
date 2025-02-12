import logging
import scanpy as sc
import os
from utils import check_matrix_files, output_adata

def pca(config):

    cache_dir = os.path.join(config['cache_dir'], 'cohort-level')
    input_path = os.path.join(cache_dir, 'normed.h5ad')
    output_path = os.path.join(cache_dir, 'pcaed.h5ad')

    if not check_matrix_files(input_path, 'pca'):
        logging.info('---------- Missing Input File: normed.h5ad ----------')        
        return

    figure_dir = os.path.join(config['figure_dir'], 'pca')
    os.makedirs(figure_dir, exist_ok=True)
    sc.settings.figdir = figure_dir
    
    adata = sc.read_h5ad(input_path)

    sc.pp.highly_variable_genes(adata, n_top_genes=2000)

    sc.tl.pca(adata)

    sc.pl.pca_variance_ratio(adata, n_pcs=50, save='pca_variance_ratio.png')

    adata.obsm['X_pca'] = adata.obsm['X_pca'][:, :config['n_pcs']]

    output_adata(adata, output_path)