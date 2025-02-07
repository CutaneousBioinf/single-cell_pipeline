import anndata as ad
from functools import partial
import logging
from multiprocessing import Pool
import os
import scanpy as sc
from utils import check_matrix_files, load_meta_data, output_adata



def process_sample(arg, input_dir):
    sample_id, harmony_by = arg[0], arg[1]

    input_path = os.path.join(input_dir, sample_id, 'removed', 'cooked_and_removed.h5ad')

    if not check_matrix_files(input_path, 'merge'):
        logging.info('---------- Missing Input File from %s ----------', sample_id)        
        return None
    
    logging.info("Merging %s...", sample_id)

    adata = sc.read_h5ad(input_path)
    adata.obs['sample_id'] = str(sample_id)
    adata.obs['harmony_by'] = str(harmony_by)

    # Updating barcodes to include sample_id
    adata.obs['barcode_sample'] = adata.obs.index
    adata.obs['barcode_sample'] = adata.obs['barcode_sample'].apply(lambda x: x + '_' + str(sample_id))
    adata.obs.index = adata.obs['barcode_sample']
    
    return adata

def parallel_Merge(config):
    meta_df = load_meta_data(config['input_dir'], config['metadata_filename'], [config['harmony_by']])
    sample_ids = meta_df.iloc[:, 0].tolist()
    batchs = meta_df.iloc[:, 1].tolist()
    args = zip(sample_ids, batchs)
    
    process_sample_partial = partial(process_sample, input_dir=config['cache_dir'])

    with Pool() as pool:
        adatas = [adata for adata in pool.map(process_sample_partial, args) if adata is not None]

    merged_adata = ad.concat(adatas)
    output_adata(merged_adata, os.path.join(config['cache_dir'], 'merged.h5ad'))