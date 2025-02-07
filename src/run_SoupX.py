from functools import partial
from utils import check_matrix_files, export_barcodes_and_genes, load_meta_data
import logging
from multiprocessing import Pool
import os
import subprocess
import warnings

warnings.filterwarnings('ignore')

def process_sample(sample_id, input_dir, output_dir):
    input_dir = os.path.join(input_dir, sample_id)

    if not check_matrix_files(input_dir, 'soupx'):
        logging.info('---------- Missing Input File from %s ----------', sample_id)        

        return
    
    logging.info("Cooking soup for %s...", sample_id)

    output_dir = os.path.join(output_dir, sample_id,  'cooked')
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate genes.tsv and barcodes.tsv
    export_barcodes_and_genes(input_dir, output_dir)

    # Generate matrix.mtx
    with open(os.devnull, 'w') as devnull:
        subprocess.run([
            "Rscript", os.path.join("src", "SoupX.R"), 
            "--input_dir", input_dir,
            "--output_dir", output_dir
        ])


def parallel_SoupX(config):
    meta_df = load_meta_data(config['input_dir'], config['metadata_filename'])
    sample_ids = meta_df.iloc[:, 0].tolist()

    process_sample_partial = partial(process_sample, input_dir=config['input_dir'], output_dir=config['cache_dir'])

    with Pool(processes=config['n_cores']) as pool:
        pool.map(process_sample_partial, sample_ids)