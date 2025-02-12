import logging
import scanpy as sc
import os
from utils import check_matrix_files, output_adata

def cluster(config):

    cache_dir = os.path.join(config['cache_dir'], 'cohort-level')
    input_path = os.path.join(cache_dir, 'umaped.h5ad')
    output_path = os.path.join(cache_dir, 'clustered.h5ad')

    sc.settings.figdir = os.path.join(config['figure_dir'], 'clustering')

    if not check_matrix_files(input_path, 'cluster'):
        logging.info('---------- Missing Input File: umaped.h5ad ----------')        
        return
    
    adata = sc.read_h5ad(input_path)

    for res in config['clustering_resolutions']:
        logging.info('Clustering with resolution: %f...', res)
        sc.tl.leiden(
            adata, 
            neighbors_key="after_harmony", 
            key_added=f"leiden_res_{res:4.2f}", 
            resolution=res, 
            flavor="igraph"
        )

    logging.info('Plotting UMAP...')
    cls = [f'leiden_res_{res:.2f}' for res in config['clustering_resolutions']]
    sc.pl.umap(
        adata,
        color=cls,
        legend_loc="on data",
        save="umap_multiple_resolution.png",
    )

    output_adata(adata, output_path)
