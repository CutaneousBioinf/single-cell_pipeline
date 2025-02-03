import anndata as ad
from functools import partial
import logging
from multiprocessing import Pool
import os
import scanpy as sc
from utils import check_matrix_files, load_meta_data, output_adata



def process_sample(arg, input_dir):
    sample_id, batch = arg[0], arg[1]

    input_path = os.path.join(input_dir, sample_id, 'removed', 'cooked_and_removed.h5ad')

    if not check_matrix_files(input_dir, 'merge'):
        logging.info('---------- Missing Input File from %s ----------', sample_id)        
        return
    
    adata = sc.read_h5ad(input_path)
    adata.obs['sample_id'] = sample_id

    # Updating barcodes to include sample_id
    adata.obs['barcode_sample'] = adata.obs.index
    adata.obs['barcode_sample'] = adata.obs['barcode_sample'].apply(lambda x: x + '_' + str(batch))
    adata.obs.index = adata.obs['barcode_sample']
    
    return adata

def parallel_Merge(config):
    meta_df = load_meta_data(config['input_dir'])
    sample_ids = meta_df.iloc[:, 0].tolist()
    batchs = meta_df.iloc[:, 1].tolist()
    args = zip(sample_ids, batchs)
    
    process_sample_partial = partial(process_sample, input_dir=config['cache_dir'])

    with Pool() as pool:
        adatas = pool.map(process_sample_partial, args)

    merged_adata = ad.concat(adatas)
    output_adata(merged_adata, os.path.join(config['cache_dir'], 'merged.h5ad'))