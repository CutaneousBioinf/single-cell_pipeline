import logging
import scanpy as sc
import os
from utils import check_matrix_files, output_adata


def normalization(config):

    cache_dir = os.path.join(config['cache_dir'], 'cohort-level')
    input_path = os.path.join(cache_dir, 'qced.h5ad')
    output_path = os.path.join(cache_dir, 'normed.h5ad')

    if not check_matrix_files(input_path, 'norm'):
        logging.info('---------- Missing Input File: qced.h5ad ----------')        
        return

    adata = sc.read_h5ad(input_path)

    sc.pp.normalize_total(adata, target_sum=1e6, inplace=True)
    sc.pp.log1p(adata, copy=False)

    output_adata(adata, output_path)