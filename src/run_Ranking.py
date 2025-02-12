import logging
import scanpy as sc
import os
from utils import check_matrix_files, output_adata

def rank(config):

    cache_dir = os.path.join(config['cache_dir'], 'cohort-level')
    input_path = os.path.join(cache_dir, 'clustered.h5ad')
    output_path = os.path.join(cache_dir, 'ranked.h5ad')

    if not check_matrix_files(input_path, 'rank'):
        logging.info('---------- Missing Input File: clustered.h5ad ----------')        
        return
    
    adata = sc.read_h5ad(input_path)

    for res in config['ranking_resolutions']:
        logging.info('Ranking with resolution: %f...', res)
        key = f"rank_genes_res_{res:.2f}"
        sc.tl.rank_genes_groups(adata, 
                                groupby=f"leiden_res_{res:.2f}", 
                                key_added=key,
                                method='wilcoxon')

        logging.info('Getting genes...')

        df = sc.get.rank_genes_groups_df(adata, key=key, group=None, pval_cutoff=config['marker_genes']['pval'], log2fc_min=config['marker_genes']['fc'])

        df.to_csv(os.path.join(cache_dir, f'marker_genes_res_{res:.2f}.csv'))
        
    output_adata(adata, output_path)

