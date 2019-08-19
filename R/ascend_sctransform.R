################################################################################
#
# ascend_sctransform.R
# description: Functions related to the normalisation of data via the
# sctransform R package
#
################################################################################

#' runSCTransform
#'
#' Normalise an \code{\linkS4class{EMSet}} with \pkg{sctransform}'s vst function
#' by Hafemeister & Satija 2019.
#'
#' @details This function is a wrapper for the vst function from the sctransform
#' package. Users may supply additional functions to the vst function.
#' 
#' @param object An \code{\linkS4class{EMSet}}
#' @param batch Normalise batches (Default: TRUE)
#' @param ngenes Number of genes to use for estimation of initial parameters 
#' (Default: 3000)
#' @param ncells Number of cells to use for the estimation of initial parameters
#' (Default: NULL)
#' @param seed Set seed (Default: 1)
#' @param verbose Print messages from vst (Default: FALSE)
#' @param ... Parameters to pass on to the vst function 
#' 
#' @return An \code{\linkS4class{EMSet}} with normalised, corrected UMI counts 
#' in the normcounts and logcounts slots, and residuals in the residuals slot.
#' 
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom SingleCellExperiment counts normcounts logcounts
#' 
#' @export
#'
runSCTransform <- function(object, batch = TRUE, seed = 1, verbose = FALSE, 
                           ngenes = 3000, ncells = NULL, ...){
  loadNamespace("sctransform")
  # Check and prepare arguments
  if (class(object) != "EMSet"){
    stop("Please supply a valid EMSet.")
  }
  
  # Check number of cells is valid
  if (!is.null(ncells) & (ncells > ncol(object))) {
    stop(sprintf("Please specify a value less than %d", ncol(object)))
  }
  
  # Prepare inputs from EMSet
  counts <- as(counts(object), "dgCMatrix")
  cell_attr <- merge(colInfo(object), colData(object), by = "cell_barcode")
  rownames(cell_attr) <- cell_attr$cell_barcode
  cell_attr <- as.data.frame(cell_attr[, c("batch", "qc_libsize", "qc_ngenes")])
  colnames(cell_attr) <- c("batch", "n_umi", "n_gene")
  
  # Run VST - load arguments into vst.args
  set.seed(seed)
  vst.args <- list(...)
  vst.args[["umi"]] <- counts
  vst.args[["batch_var"]] <- "batch"
  vst.args[["cell_attr"]] <- cell_attr
  vst.args[["show_progress"]] <- verbose
  vst.args[['return_cell_attr']] <- TRUE
  vst.args[['return_gene_attr']] <- TRUE
  vst.args[['return_corrected_umi']] <- TRUE
  vst.args[['n_cells']] <- ncells
  vst.args[["n_genes"]] <- ngenes
  
  # Run VST
  vst_outs <- do.call(what = "vst", args = vst.args)
  
  # Extract values from VST
  gene_attr <- vst_outs$gene_attr
  cell_attr <- vst_outs$cell_attr
  
  # Sort variance
  var_df <- data.frame(gene_id = rownames(gene_attr), 
                       residual_variance = gene_attr[, "residual_variance"])
  var_df <- var_df[order(var_df$residual_variance, decreasing = TRUE), ]
  var_df$ranking <- 1:nrow(var_df)
  var_df <- var_df[match(rownames(gene_attr), var_df$gene_id), ]
  gene_attr$var_ranking <- var_df$ranking
  
  # Extract matrices
  norm_counts <- vst_outs$umi_corrected
  logcounts_matrix <- log2(norm_counts + 1)
  res_matrix <- vst_outs$y
  
  # Regress controls
  object <- object[rownames(norm_counts), colnames(norm_counts)]
  normcounts(object) <- as(norm_counts, "dgCMatrix")
  logcounts(object) <- as(logcounts_matrix, "dgCMatrix")
  regcounts(object) <- as.matrix(res_matrix)
  
  row_info <- rowInfo(object)
  row_data <- rowData(object)
  col_info <- colInfo(object)
  col_data <- colData(object)
  
  row_data <- cbind(row_data, S4Vectors::DataFrame(gene_attr))
  col_data <- cbind(col_data, S4Vectors::DataFrame(cell_attr))
  
  colData(object) <- col_data
  rowData(object) <- row_data
  
  log <- progressLog(object)
  log$NormalisationMethod <- "SCT"
  log$SCTransform <- vst.args[which(!(names(vst.args) %in% c("umi", "cell_attr")))]
  progressLog(object) <- log
  return(object)
}
  