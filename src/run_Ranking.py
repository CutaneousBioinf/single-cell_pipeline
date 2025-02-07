import logging
import scanpy as sc
import os
from utils import check_matrix_files, output_adata

def rank(config):

    cache_dir = config['cache_dir']
    input_path = os.path.join(cache_dir, 'clustered.h5ad')
    output_path = os.path.join(cache_dir, 'ranked.h5ad')

    if not check_matrix_files(input_path, 'rank'):
        logging.info('---------- Missing Input File: normed.h5ad ----------')        
        return
    
    adata = sc.read_h5ad(input_path)

    for res in config['resolutions']:
        logging.info('Ranking with resolution: %f...', res)
        sc.tl.rank_genes_groups(adata, 
                                groupby=f"leiden_res_{res:.2f}", 
                                key_added=f"rank_genes_res_{res:.2f}")
        
    output_adata(adata, output_path)

