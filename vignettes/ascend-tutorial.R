## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load_expression-----------------------------------------------------
library(ascend)
data(raw_set)
EMSet <- raw_set
counts <- counts(EMSet)
counts[1:5, 1:5]

## ----cell_info_knitr-----------------------------------------------------
# Load dataframe stored in colInfo slot.
col_info <- colInfo(EMSet)
print(col_info)

## ----addTHY1-------------------------------------------------------------
# Retrieve list of batches stored in colInfo
batch_list <- col_info$batch

# Replace 1 values in col_info with "Positive" by identifying indices of 1 
# values in the list
batch_list[which(batch_list == 1)] <- "Positive"

# Replace 2 values in col_info with "Negative" by identifying indices of 2 
# values in the list
batch_list[which(batch_list == 2)] <- "Negative"

# Add as a new column in col_info data frame.
col_info$THY1 <- batch_list

# Update EMSet with it
colInfo(EMSet) <- col_info

# Print colInfo 
print(colInfo(EMSet))

## ----labelCells----------------------------------------------------------
gene_markers <- c("POU4F1", "POU4F2", "POU4F3")
EMSet <- addGeneLabel(EMSet, gene = gene_markers)
colInfo(EMSet)

## ----gene_info_knitr-----------------------------------------------------
# Load dataframe stored in colInfo slot.
row_info <- rowInfo(EMSet)
print(row_info)

## ----switchAnnotation----------------------------------------------------
ensembl_set <- convertGeneID(EMSet, new.annotation = "ensembl_gene_id")
print("rowInfo dataframe")
rowInfo(ensembl_set)

print("Count matrix")
counts(ensembl_set)[1:5, 1:5]

## ----set_controls--------------------------------------------------------
# Example controls - generate them from row names
controls <- list(Mt = grep("^Mt-", rownames(EMSet), ignore.case = TRUE, value = TRUE),
                 Rb = grep("^Rps|^Rpl", rownames(EMSet), ignore.case = TRUE, value = TRUE))

# Set controls in the EMSet                 
controls(EMSet) <- controls

# Retrieve controls
print(controls(EMSet))

## ----remove_controls, eval = FALSE---------------------------------------
#  EMSet <- excludeControl(EMSet, control = c("Mt", "Rb"))

## ----load_data, eval = FALSE---------------------------------------------
#  # Get sample data to make a new set
#  counts <- counts(EMSet)
#  col_info <- colInfo(EMSet)
#  row_info <- rowInfo(EMSet)
#  controls <- controls(EMSet)
#  
#  # Creating an EMSet from scratch
#  em_set <- EMSet(list(counts = counts),
#                     colInfo = col_info,
#                     rowInfo = row_info,
#                     controls = controls)
#  
#  # Loading an EMSet from a SingleCellExperiment
#  single_cell_experiment <- SingleCellExperiment(list(counts = counts))
#  em_set <- EMSet(single_cell_experiment,
#                     colInfo = col_info,
#                     rowInfo = row_info,
#                     controls = controls)
#  

## ----get_set_examples----------------------------------------------------
# Load a 'complete' EMSet
data(analyzed_set)
EMSet <- analyzed_set

# Get assays
count_matrix <- counts(EMSet)
norm_matrix <- normcounts(EMSet)
logcounts_matrix <- logcounts(EMSet)

# Set assays
counts(EMSet) <- count_matrix
normcounts(EMSet) <- norm_matrix
logcounts(EMSet) <- logcounts_matrix

# Get gene and cell information
col_info <- colInfo(EMSet)
col_data <- colData(EMSet)
row_info <- rowInfo(EMSet)
row_data <- rowData(EMSet)

# Set gene and cell information
colInfo(EMSet) <- col_info
colData(EMSet) <- col_data
rowInfo(EMSet) <- row_info
rowData(EMSet) <- row_data

# Get reduced dimensionality data
tsne_matrix <- reducedDim(EMSet, "TSNE")
pca_matrix <- reducedDim(EMSet, "PCA")
umap_matrix <- reducedDim(EMSet, "UMAP")

# Get cluster analysis
clusterAnalysis <- clusterAnalysis(EMSet)

# Get progress log
progressLog <- progressLog(EMSet)

## ----dataframe_accesors--------------------------------------------------
# Reduce EMSet to first ten cells and genes
tiny_EMSet <- EMSet[1:10,1:10]

# Review content in smaller dataset 
print(counts(tiny_EMSet))
print(colInfo(tiny_EMSet))
print(rowInfo(tiny_EMSet))

## ----subsetCondition-----------------------------------------------------
# Subset batch 1 from the combined EMSet
Batch1_EMSet <- subsetCondition(EMSet, by = "batch", conditions = list(batch = c(1)))

## ----update_object, eval = FALSE-----------------------------------------
#  # Read in old EMSet stored as an RDS file
#  legacy_EMSet <- readRDS("LegacyEMSet.rds")
#  
#  # Update to new object
#  # Please ensure your new object has the same name as the old object
#  legacy_EMSet <- updateObject(legacy_EMSet)
#  

## ----runTSNE-------------------------------------------------------------
EMSet <- runTSNE(EMSet, dims = 2, PCA = TRUE, seed = 1)

## ----reducedDim----------------------------------------------------------
# Retrieve dimensionality-reduced matrices
tsne_matrix <- reducedDim(EMSet, "TSNE")
pca_matrix <- reducedDim(EMSet, "PCA")
umap_matrix <- reducedDim(EMSet, "UMAP")

# Display retrieved matrices
print("t-SNE matrix")
print(tsne_matrix[1:10, ])

print("PCA matrix")
print(pca_matrix[1:10, 1:5])

print("UMAP matrix")
print(umap_matrix[1:10, ])


## ----plotTSNE, fig.width = 6, fig.height = 5-----------------------------
tsne_plot <- plotTSNE(EMSet, group = "cluster")
tsne_plot

## ----colour_tsne, fig.width = 6, fig.height = 5--------------------------
library(ggplot2)
# Set colour palette
tsne_plot <- tsne_plot + scale_color_manual(values=c("#bb5f4c", 
                                                     "#8e5db0", 
                                                     "#729b57"))

# Add title
tsne_plot <- tsne_plot + ggtitle("Clusters", subtitle = "tSNE plot") 

# Put legend on the bottom
tsne_plot <- tsne_plot + theme(legend.position = "bottom")

tsne_plot

## ----umap_plot, fig.width = 6, fig.height = 5----------------------------
library(ggplot2)

# Generate plot using plotUMAP
umap_plot <- plotUMAP(EMSet, group = "cluster", Dim1 = 1, Dim2 = 2)

# Set colour palette
umap_plot <- umap_plot + scale_color_manual(values=c("#bb5f4c", 
                                                     "#8e5db0", 
                                                     "#729b57"))

# Add title
umap_plot <- umap_plot + ggtitle("Clusters", subtitle = "UMAP plot") 

# Put legend on the bottom
umap_plot <- umap_plot + theme(legend.position = "bottom")
umap_plot

## ----convert2SCE---------------------------------------------------------
# Convert to SingleCellExperiment
SingleCellExperiment <- convert(EMSet, to = "sce")

head(SingleCellExperiment)

# Revert to SingleCellExperiment
EMSet <- convert(SingleCellExperiment, to = "EMSet")
head(EMSet)


## ----convertToSeurat, eval = FALSE---------------------------------------
#  SeuratObject <- convert(EMSet, to = "seurat")
#  SeuratObject
#  EMSet <- convert(SeuratObject, to = "EMSet")

## ----convertToScone, eval = FALSE----------------------------------------
#  sconeObject <- convert(EMSet, to = "scone")
#  sconeObject
#  EMSet <- convert(sconeObject, to = "EMSet")

## ----scran_normalisation, eval = FALSE-----------------------------------
#  scran_normalised <- scranNormalise(EMSet, quickCluster = FALSE, min.mean = 1e-05)

## ----scran_cellcycle, eval = FALSE---------------------------------------
#  # Convert identifiers to ENSEMBL gene identifiers
#  ensembl_set <- convertGeneID(RawSet, new.annotation = "ensembl_gene_id")
#  
#  # Load training data from scran
#  training_data <- readRDS(system.file("exdata", "human_cycle_markers.rds",
#                                       package = "scran"))
#  
#  # Run cell cycle
#  scran_cellcycle <- scranCellCycle(ensembl_set, training_set = training_data)
#  
#  # Show cell cycle results
#  colInfo(scran_cellcycle)

## ----DESeq, eval = FALSE-------------------------------------------------
#  cluster2_vs_others <- runDESeq(EMSet, group = "cluster",
#                                 condition.a = 2, condition.b = c(1, 3),
#                                 ngenes = 1500, parallel = TRUE)

## ----DESeq2, eval = FALSE------------------------------------------------
#  cluster1_vs_others <- runDESeq2(EMSet, group = "cluster",
#                                 condition.a = 1, condition.b = c(2, 3),
#                                 ngenes = 1500)

