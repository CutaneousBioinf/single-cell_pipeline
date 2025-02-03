import logging
import scanpy as sc
import os
from utils import check_matrix_files, output_adata

def clustering(config):

    cache_dir = config['cache_dir']
    input_path = os.path.join(cache_dir, 'harmonyed.h5ad')
    output_path = os.path.join(cache_dir, 'clustered.h5ad')

    sc.settings.figdir = config['figure_dir']

    if not check_matrix_files(input_path, 'cluster'):
        logging.info('---------- Missing Input File: normed.h5ad ----------')        
        return
    
    adata = sc.read_h5ad(input_path)

    logging.info('Clustering with before harmony data')
    sc.pp.neighbors(adata, key_added='before_harmony', use_rep='X_pca')
    sc.tl.umap(adata, neighbors_key='before_harmony')
    adata.obsm['X_umap_before_harmony'] = adata.obsm['X_umap'].copy()
    sc.pl.umap(adata, neighbors_key='before_harmony', color=config['harmony_by'], title='UMAP before Harmony', save='umap_before_harmony.png')

    logging.info('Clustering with after harmony data')
    sc.pp.neighbors(adata, key_added='after_harmony', use_rep='X_pca_harmony')
    sc.tl.umap(adata, neighbors_key='after_harmony')
    adata.obsm['X_umap_after_harmony'] = adata.obsm['X_umap'].copy()
    sc.pl.umap(adata, neighbors_key='after_harmony', color=config['harmony_by'], title='UMAP after Harmony', save='umap_after_harmony.png')  

    output_adata(adata, output_path)
