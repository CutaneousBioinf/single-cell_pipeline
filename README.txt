Welcome to the Single-Cell Pipeline!


# 0. Pipeline Overview

    Check out the second page of the slides in the repo


# 1. Setting up your local repo:

    $ git clone https://github.com/lica-2025/single-cell_pipeline.git
    $ cd single-cell_pipeline


# 2. Setting up environment:

    $ mamba env create -f environment.yml
    $ mamba activate sc-pipeline
    (you can change the environment name in the first line of environment.yml)


# 3. Take a look at the repo directory:

    path_to_your_directory/sc-pipeline/
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
        │   └── run_Clustering.py
        ├── config/
        │   ├── config_sample.json
        │   └── (config.json)
        ├── log/
        │   ├── (info.log)
        │   ├── (err.log)
        │   └── (pipeline_yyyymmdd_hhmmss.pid)
        ├── README.txt
        └── run_sc-pipeline.sh


# 4. Setting up the configuration file:

    Make a copy of the sample configuration file
        $ cp /path_to_your_directory/sc-pipeline/config/config_sample.json /path_to_your_directory/sc-pipeline/config/config.json 

    {

        "input_dir": path_to_input_directory/,
            ---------- please read the next bullet point for details ----------

        "cache_dir": path_to_cache_directory/,
            directory where temporary output objects from each step are saved
            no need to pre-create

        "figure_dir": path_to_figure_directory/,
            directory where figures from the pipelien are saved
            no need to pre-create

        "harmony_by":
            should be one of the columns names in adata.obs
            this column will be used in batch-correction

        "n_cores": 
            #cores used for the pipeline 
            this will only affect cooking soup, removing doublet, and merging steps
            as these are running are sample level thus being parallelized
    }


# 5. Preparing input files

    /path_to_input_directory/
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


# 6. Running sc-pipeline

        $ bash run_sc-pipeline.sh


# 7. Monitoring

    Is the job still running?
        $ ps -u (your_server_login_name)
        check if the job's pid is still there

    How to find the job's pid?
        check /log/pipeline_yyyymmdd_hhmmss.pid, where the time in the filename is the starting time (GMT)

    Which step is the job running now?
        check /log/info.log for current progress

    What's wrong with the job?
        check /log/err.log for any interuptions

    How to terminate the pipeline?
        $ kill -9 (pid)
        there might still be sub-threads running, to kill them:
        $ pkill -u (your_server_login_name) python
        $ pkill -u (your_server_login_name) R
        Cautious: these will kill all your other python and R jobs (if exists)


# 8. Output

    find figures in /path_to_figure_directory/

    find objects in /path_to_cache_directory/

    the final object is at /path_to_cache_directory/clustered.h5ad
