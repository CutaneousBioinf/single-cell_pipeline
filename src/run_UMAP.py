import logging
import scanpy as sc
import os
from utils import check_matrix_files, output_adata

def umap(config):

    cache_dir = os.path.join(config['cache_dir'], 'cohort-level')
    input_path = os.path.join(cache_dir, 'harmonyed.h5ad')
    output_path = os.path.join(cache_dir, 'umaped.h5ad')

    sc.settings.figdir = os.path.join(config['figure_dir'], 'harmony')

    if not check_matrix_files(input_path, 'umap'):
        logging.info('---------- Missing Input File: harmonyed.h5ad ----------')        
        return
    
    adata = sc.read_h5ad(input_path)

    logging.info('Running UMAP with before harmony data')
    sc.pp.neighbors(adata, key_added='before_harmony', use_rep='X_pca')
    sc.tl.umap(adata, neighbors_key='before_harmony')
    adata.obsm['X_umap_before_harmony'] = adata.obsm['X_umap'].copy()
    sc.pl.umap(adata, neighbors_key='before_harmony', color='harmony_by', title='UMAP before Harmony', save='umap_before_harmony.png')

    logging.info('Running UMAP with after harmony data')
    sc.pp.neighbors(adata, key_added='after_harmony', use_rep='X_pca_harmony')
    sc.tl.umap(adata, neighbors_key='after_harmony')
    # adata.obsm['X_umap_after_harmony'] = adata.obsm['X_umap'].copy()
    sc.pl.umap(adata, neighbors_key='after_harmony', color='harmony_by', title='UMAP after Harmony', save='umap_after_harmony.png')  

    output_adata(adata, output_path)
