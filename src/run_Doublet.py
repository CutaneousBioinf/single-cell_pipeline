from functools import partial
import logging
from multiprocessing import Pool
import os
import scanpy as sc
from utils import check_matrix_files, load_meta_data, read_10x_anndata, output_adata

def process_sample(sample_id, input_dir, output_dir):
    input_dir = os.path.join(input_dir, sample_id, 'cooked')

    if not check_matrix_files(input_dir, type='doublet'):
        logging.info('---------- Missing Input File from %s ----------', sample_id)        

        return

    logging.info('Removing doublet for %s', sample_id)

    output_dir = os.path.join(output_dir, sample_id, 'removed')
    os.makedirs(output_dir, exist_ok=True)

    adata = read_10x_anndata(input_dir)

    # add doublet info to adata.obs, not actually removing doublets
    sc.pp.scrublet(adata)
    
    output_adata(adata, os.path.join(output_dir,  'cooked_and_removed.h5ad'))

def parallel_Doublet(config):
    meta_df = load_meta_data(config['input_dir'], config['metadata_filename'])
    sample_ids = meta_df.iloc[:, 0].tolist()

    process_sample_partial = partial(process_sample, input_dir=config['cache_dir'], output_dir=config['cache_dir'])

    with Pool(processes=config['n_cores']) as pool:
        pool.map(process_sample_partial, sample_ids)
