import os 
import pandas as pd
import scanpy as sc

def check_matrix_files(input_dir, type):
    if type == 'soupx':
        filt_matrix_path = os.path.join(input_dir, 'filtered_feature_bc_matrix.h5')
        raw_matrix_path = os.path.join(input_dir, 'raw_feature_bc_matrix.h5')
        return os.path.exists(filt_matrix_path) and os.path.exists(raw_matrix_path)
    elif type == 'doublet':
        count_path = os.path.join(input_dir, 'matrix.mtx')
        genes_path = os.path.join(input_dir, 'genes.tsv')
        barcodes_path = os.path.join(input_dir, 'barcodes.tsv')
        return os.path.exists(count_path) and os.path.exists(genes_path) and os.path.exists(barcodes_path)
    else:
        return os.path.exists(input_dir)


def export_barcodes_and_genes(input_dir, output_dir):
    file_path = os.path.join(input_dir, 'filtered_feature_bc_matrix.h5')
    adata = sc.read_10x_h5(file_path, gex_only=True)
    genes = pd.DataFrame({'gene_names': adata.var.index}, index=adata.var['gene_ids'])
    genes.to_csv(os.path.join(output_dir, 'genes.tsv'), header=False, sep='\t')
    adata.obs[[]].to_csv(os.path.join(output_dir, 'barcodes.tsv'), header=False, sep='\t')

def load_meta_data(input_dir, filename='metadata.csv', extra_columns=[]):
    # In run_Merge.py, if adds more columns to extra_columns, 
    # make sure config['harmony_by'] is the first one.
    meta_path = os.path.join(input_dir, filename)
    meta_df = pd.read_csv(meta_path)
    meta_df = meta_df[['Sample_ID'] + extra_columns]
    return meta_df

def output_adata(adata, output_path):
    adata.write_h5ad(filename=output_path)

def read_10x_anndata(file_path):
    adata = sc.read_10x_mtx(file_path, var_names='gene_symbols', cache=False)
    return adata