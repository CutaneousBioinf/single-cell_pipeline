# Welcome to the Single-Cell Pipeline!


## 0. Pipeline Overview

Check out the second page of the slides in the repo


## 1. Setting up your local repo:
    $ git clone https://github.com/lica-2025/single-cell_pipeline.git
    $ cd single-cell_pipeline


## 2. Setting up environment:
    $ mamba env create -f environment.yml
    $ mamba activate sc-pipeline

You can change the environment name in the first line of environment.yml


## 3. Take a look at the repo directory:

    path_to_your_directory/single-cell_pipeline/
        ├── src/
        │   ├── main.py
        │   ├── utils.py
        │   ├── run_SoupX.py
        │   ├── SoupX.R
        │   ├── run_Doublet.py
        │   ├── run_Merge.py
        │   ├── run_QC.py
        │   ├── run_Norm.py
        │   ├── run_PCA.py
        │   ├── run_Harmony.py
        │   └── run_UMAP.py
        │   └── run_Clustering.py
        │   └── run_Ranking.py
        │   └── run_Annotation.py
        ├── config/
        │   ├── config_sample.json
        │   └── (config.json)
        ├── (log/)
        │   ├── (info.log)
        │   ├── (err.log)
        │   └── (pipeline_yyyymmdd_hhmmss.pid)
        ├── run_sc-pipeline.sh
        ├── environment.yml
        ├── .gitignore
        ├── README.md
        └── LICENSE



## 4. Setting up the configuration file:

Make a copy of the sample configuration file

    $ cp /path_to_your_directory/sc-pipeline/config/config_sample.json /path_to_your_directory/sc-pipeline/config/config.json 


### Directories

#### input_dir
Path to the input directory, please read the next bullet point for details.

#### cache_dir
Path to the cache directory, where temporary output objects and csvs from each step are saved. No need to pre-create.

#### figure_dir
Path to the figure directory, where figures from the pipeline are saved. No need to pre-create.

#### metadata_filename
You can change this to match the metadata file's name in the input directory.

### Parameters

#### qc_metrics
Dictionary dictating QC metrics. Supported keys: min_gene_per_cell, max_doublet_score, max_mt_percentage.

#### n_pcs
Number of Principal Components to keep at the end of the PCA step.

#### harmony_by
Should be one of the columns names in the metadata. The column will be used in batch-correction.

#### marker_genes
Dictionary dictating marker gene identifying metrics. Should include only: pval (max p-value) and fc (min fold change). fc should be positive number. Only up-regulated genes will be identified.

#### annotation
Dictionary assigning cell type to each cluster. Keep empty if don't want to run annotation, otherwise make sure length of this dictionary matches #cluster with resolution == annotation_resolution.

### Resolutions

#### clustering_resolutions
Plot UMAP with different resolutions to decide the best resolution. Greater value of resolution generates more clusteres.

#### ranking_resolutions
Identify significant genes with different resolutions. ranking_resolutions must be a subset of clustering_resolutions because ranking is based on clustering.

##### annotation_resolution
Plot UMAP with cluster names. annotation_resolutions must be an element from ranking_resolutions because annotation is based on clustering.

### Performance

#### n_cores
This is #cores used for the pipeline, which will only affect cooking soup, removing doublet, and merging steps, as these are running on sample level thus being parallelized.



## 5. Preparing input files

    path_to_input_directory/
    ├── metadata.csv
    │       Sample_ID,Batch
    │       sample_1,1
    │       sample_2,1
    │       sample_3,2
    │       sample_4,3
    │       sample_5,3
    │       sample_6,3
    │       ...
    ├── sample_1/
    │       ├── raw_feature_bc_matrix.h5
    │       └── filtered_feature_bc_matrix.h5
    ├── sample_2/
    │       ├── raw_feature_bc_matrix.h5
    │       └── filtered_feature_bc_matrix.h5
    ├── .../


## 6. Running sc-pipeline

    $ bash run_sc-pipeline.sh


## 7. Monitoring

### Is the job still running?

    $ ps -u (your_server_login_name)

### How to find the job's pid?
Check 

    /log/pipeline_yyyymmdd_hhmmss.pid

where the time in the filename is the starting time (GMT)

### Which step is the job running now?
Check 

    /log/info.log
    
for current progress

### What's wrong with the job?
Check 

    /log/err.log
    
for any interuptions

### How to terminate the pipeline?

    $ kill -9 (pid)

There might still be sub-threads running, to kill them:

    $ pkill -u (your_server_login_name) python
    $ pkill -u (your_server_login_name) R

Cautious: these will kill all your other python and R jobs (if exists)


## 8. Output

Find figures in 

    path_to_figure_directory/

Find objects in

    path_to_cache_directory/

The final object is at

    path_to_cache_directory/clustered.h5ad
