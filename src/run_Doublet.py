from functools import partial
import logging
from multiprocessing import Pool
import os
import scanpy as sc
from utils import check_matrix_files, load_meta_data, output_adata, read_10x_anndata

def process_sample(sample_id, cache_dir):
    cache_dir = os.path.join(cache_dir, 'sample-level', sample_id)

    if not check_matrix_files(cache_dir, type='doublet'):
        logging.info('---------- Missing Input File from %s ----------', sample_id)        

        return

    logging.info('Removing doublet for %s', sample_id)

    adata = read_10x_anndata(cache_dir)

    # add doublet info to adata.obs, not actually removing doublets
    sc.pp.scrublet(adata)
    
    output_adata(adata, os.path.join(cache_dir, "doublet_calculated.h5ad"))

def parallel_Doublet(config):
    meta_df = load_meta_data(config['input_dir'], config['metadata_filename'])
    sample_ids = meta_df.iloc[:, 0].tolist()

    process_sample_partial = partial(process_sample, cache_dir=config['cache_dir'])

    with Pool(processes=config['n_cores']) as pool:
        pool.map(process_sample_partial, sample_ids)
