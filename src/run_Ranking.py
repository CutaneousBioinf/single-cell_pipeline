import logging
import scanpy as sc
import os
from utils import check_matrix_files, output_adata

def rank(config):

    cache_dir = config['cache_dir']
    input_path = os.path.join(cache_dir, 'clustered.h5ad')
    output_path = os.path.join(cache_dir, 'ranked.h5ad')

    if not check_matrix_files(input_path, 'rank'):
<<<<<<< HEAD
        logging.info('---------- Missing Input File: clustered.h5ad ----------')        
=======
        logging.info('---------- Missing Input File: normed.h5ad ----------')        
>>>>>>> 2bbac3071cc55c1d7853a1d1a92e810372a17376
        return
    
    adata = sc.read_h5ad(input_path)

    for res in config['resolutions']:
        logging.info('Ranking with resolution: %f...', res)
        sc.tl.rank_genes_groups(adata, 
                                groupby=f"leiden_res_{res:.2f}", 
<<<<<<< HEAD
                                key_added=f"rank_genes_res_{res:.2f}",
                                method='wilcoxon')

    logging.info('Ranking with resolution: %f...', config['resolutions'][-1])
    res = f"rank_genes_res_{config['resolutions'][-1]:.2f}"

    df = sc.get.rank_genes_groups_df(adata, key=res, group=None, pval_cutoff=0.05, log2fc_min=1.0)

    df.to_csv(os.path.join(cache_dir, 'markers.csv'))
=======
                                key_added=f"rank_genes_res_{res:.2f}")
>>>>>>> 2bbac3071cc55c1d7853a1d1a92e810372a17376
        
    output_adata(adata, output_path)

