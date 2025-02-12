import logging
import scanpy as sc
import os
from utils import check_matrix_files, output_adata


def harmony(config):

    cache_dir = os.path.join(config['cache_dir'], 'cohort-level')
    input_path = os.path.join(cache_dir, 'pcaed.h5ad')
    output_path = os.path.join(cache_dir, 'harmonyed.h5ad')

    if not check_matrix_files(input_path, 'harmony'):
        logging.info('---------- Missing Input File: pcaed.h5ad ----------')        
        return
    
    adata = sc.read_h5ad(input_path)

    sc.external.pp.harmony_integrate(adata, key='harmony_by')

    output_adata(adata, output_path)