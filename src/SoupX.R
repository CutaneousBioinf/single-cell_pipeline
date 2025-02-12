library(Seurat)
library(Matrix)
library(SoupX)
library(argparse)

parser <- ArgumentParser()

parser$add_argument("--input_dir", type = "character", required = TRUE, 
                    help = "Path to the input directory")
parser$add_argument("--output_dir", type = "character", required = TRUE, 
                    help = "Path to the output directory")

args <- parser$parse_args()
input_dir <- args$input_dir
output_dir <- args$output_dir

filt_matrix_path <- file.path(input_dir, 'filtered_feature_bc_matrix.h5')
raw_matrix_path <- file.path(input_dir, 'raw_feature_bc_matrix.h5')

raw.matrix <- Read10X_h5(raw_matrix_path)
filt.matrix <- Read10X_h5(filt_matrix_path)

seurat <- CreateSeuratObject(counts = filt.matrix)


# ---------- You might want to remove empty genes at this step ----------
# ---------- Redundant features will slowdown (if not disable) the pipeline ----------

# --- sometimes input file has multiple layers, which will also include extra features ---
# --- In the following example, we only use the couts.Gene layer and the extra featrues are automatically removed ---
filt.matrix <- LayerData(seurat[["RNA"]], layer = "counts.Gene Expression") # get the layer we want
seurat <- CreateSeuratObject(counts = filt.matrix)  # re-create the object using one-layer matrix
seurat_raw <- CreateSeuratObject(counts = raw.matrix)   # same steps for raw matrix
raw.matrix <- LayerData(seurat_raw[["RNA"]], layer = "counts.Gene Expression")

seurat <- NormalizeData(seurat, verbose = F)
seurat <- FindVariableFeatures(seurat, verbose = F)
seurat <- ScaleData(seurat, verbose = F)
seurat <- RunPCA(seurat, verbose = F)
seurat <- FindNeighbors(seurat, verbose = F)
seurat <- FindClusters(seurat, verbose = F)

# saveRDS(seurat, file = file.path(output_dir, "seurat_clustered_not_soup.rds"))

soup.channel  <- SoupChannel(raw.matrix, filt.matrix)
soup.channel  <- setClusters(soup.channel, setNames(seurat@meta.data$seurat_clusters, rownames(seurat@meta.data)))
soup.channel  <- autoEstCont(soup.channel)
adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)

writeMM(adj.matrix, file.path(output_dir, "matrix.mtx"))