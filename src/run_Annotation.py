import logging
import scanpy as sc
import os
from utils import check_matrix_files, output_adata

def annotation(config):

    cache_dir = os.path.join(config['cache_dir'], 'cohort-level')
    input_path = os.path.join(cache_dir, 'ranked.h5ad')
    output_path = os.path.join(cache_dir, 'annotated.h5ad')

    sc.settings.figdir = os.path.join(config['figure_dir'], 'annotation')

    if not check_matrix_files(input_path, 'rank'):
        logging.info('---------- Missing Input File: ranked.h5ad ----------')        
        return
    
    adata = sc.read_h5ad(input_path)

    res = config['annotation_resolution']
    logging.info('Annotating with resolution: %f...', res)

    annotation_key = f'annotation_res_{res:.2f}'
    adata.obs[annotation_key] = adata.obs[f"leiden_res_{res:.2f}"].map(config['annotation'])

    sc.pl.umap(
            adata,
            color=annotation_key,
            legend_loc="on data",
            title=f"Cluster Annotation, res={res:.2f}",
            save = annotation_key + '.png'
    )

    output_adata(adata, output_path)

