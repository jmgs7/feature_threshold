#' @title .ComputeFeatureThresholdSeurat
#'
#' @description
#' Computes an automatic threshold on normalized feature expression to separate low-quality
#' droplets or empty droplets from real cells. This function is intented to follow the
#' MALAT1 threshold approach developed by BaderLab (see references). This function is
#' designed to work with Seurat objects and assumes that the input Seurat object has been
#' normalized and contains a the searched feature in its data layers.
#'
#' @details
#' MALAT1 is a nuclear-retained lncRNA whose expression correlates strongly with
#' nuclear (intronic) RNA content and splice ratio, making it a convenient proxy
#' for cell quality in droplet-based scRNA-seq.
#' Low MALAT1 values typically correspond to empty droplets, cytosolic debris,
#' or damaged cells, whereas higher values correspond to intact nuclei.
#' This function assumes you provide a SeuratObject with normalized counts
#' and a "MALAT1" feature in its data layers. It computes a threshold to
#' flag low-quality cells based on the distribution of MALAT1 expression.
#' The computation is performed per provided layer (e.g., "data.pool1", "data.pool2", etc.).
#' The threshold deterined per layer and a boolean value indicating whether each cell passes
#' the threshold are stored in the Seurat object's metadata.
#'
#' @param SeuratObject A SeuratObject with normalized counts and a the searched feature in its data layers.
#' @param assay Character string specifying the assay in the SeuratObject to use for MALAT1 expression. Default is "RNA".
#' @param layers Character vector specifying the layers in the SeuratObject to process.
#'   Default is `NULL`, which processes all layers containing "data" in their names.
#' @param feature Character string specifying the feature to use for MALAT1 expression. Default is "MALAT1". This
#'   allows to search for MALAT1 in case is encoded with other symbol and to use the function for other features.
#' @param ... Additional arguments passed to the underlying `.ComputeFeatureThreshold` function, such as:
#'  \itemize{
#'    \item \code{bw.bandwidth}: Bandwidth for kernel density estimation (default: 0.01).
#'    \item \code{chosen.min}: Chosen minimum which a peak should be considered the dataset peak (default: 2).
#'    \item \code{smooth.spar}: Smoothing parameter for density estimation (default: 2).
#'    \item \code{abs.min}: Absolute minimum threshold (default: 1).
#'    \item \code{rough.max}: Rough expected position of the MALAT1 expression peak (default: 6).
#'    \item \code{conservative.threshold}: Conservative threshold to apply when impossible to find local minimum (default: 2).
#' }
#'
#' @import Seurat
#' @import SeuratObject
#' @importFrom data.table rbindlist
#'
#'
#' @return:
#' \describe{
#'   \item{SeuratObject} {The input SeuratObject with the MALAT1 threshold applied.}
#' }
#'
#'@noRd
#'
#' @references
#' Clarke, Bader et al. "MALAT1 expression indicates cell quality in
#' single-cell RNA sequencing data." bioRxiv (2024).
#' \url{https://www.biorxiv.org/content/10.1101/2024.07.14.603469v2}
#'
#' @examples
#' \dontrun{
#'   threshold.res <- .CalculateFeatureThresholdSeurat(
#'     SeuratObject = seurat_obj,
#'     Assay = "RNA",
#'     Layers = c("data.pool1", "data.pool2"),
#'     feature = "MALAT1",
#'     bw.bandwidth = 0.01,
#'     chosen.min = 2,
#'     smooth.spar = 2,
#'     abs.min = 1,
#'     rough.max = 6,
#'     conservative.threshold = 2
#'   )
#' }
#'
#' @seealso
#' \itemize{
#'   \item DropletQC: nuclear fraction-based QC for empty droplets and damaged cells.
#'   \item EmptyDrops: ambient RNA-based empty droplet detection.
#' }

.CalculateFeatureThresholdSeurat <- function(
  SeuratObject,
  assay = "RNA",
  layers = NULL,
  feature = "MALAT1",
  ...
) {
  # Check if the specified assay exists in the Seurat object
  if (!assay %in% names(SeuratObject@assays)) {
    stop(paste("Assay", assay, "not found in the Seurat object."))
  }

  # Fetch all layers present in the Seurat object for the specified assay
  object.layers <- SeuratObject::Layers(
    SeuratObject,
    assay = assay,
  )

  # Check if the user provided layers are present in the Seurat object.
  if (!is.null(layers)) {
    if ((!all(layers %in% object.layers))) {
      stop(
        "Some specified layers are not present in the Seurat object. Please check the provided layers argument."
      )
    }
  } else {
    # If no layers are provided, check if any log-normalized data layers are present in the Seurat object.
    if (!any(grepl("data", object.layers))) {
      stop(
        "User did not specified any layers and no log-normalized data layers found in the Seurat object. Please specify layers to process."
      )
    }
    # If no layers are provided, fetch all log-normalized 'data' layers names from the Seurat object.
    layers <- object.layers[grepl("data", object.layers)]
  }

  # Check if the Seurat object contains a "MALAT1" feature in the specified assay.
  # If not, stop the function and return an error message.
  tryCatch(
    SeuratObject::FetchData(SeuratObject, assay = assay, vars = feature),
    error = function(e) {
      stop(
        "The Seurat object does not contain a '",
        feature,
        "' feature in the specified assay. Please ensure that the '",
        feature,
        "' feature is present. Error details from Seurat::FetchData: ",
        e$message
      )
    }
  )

  # The results of the feature threshold computation will be stored in the Seurat object's metadata.
  feature.metadata <- lapply(
    layers,
    function(layer) {
      # Retrieve the log-normalized MALAT1 expression values for the current layer.
      # We use suppressWarnings to avoid the warnings that arise when fetching data f
      # from specific layers.
      feature.logcounts <- suppressWarnings(SeuratObject::FetchData(
        SeuratObject,
        assay = assay,
        vars = feature,
        layer = layer
      ))

      # Compute the MALAT1 threshold using the .ComputeFeatureThreshold function, passing any additional arguments.
      threshold.value <- round(
        .ComputeFeatureThreshold(feature.logcounts[[feature]], ...),
        digits = 2
      )

      # We will store the results in a data frame with cell IDs, computed threshold,
      # and whether each cell passes the threshold.
      # The cells id are the row names of the feature.logcounts data frame to ensure
      # proper alignment with the Seurat object's metadata.
      cell.id <- row.names(feature.logcounts)
      # Create a vector of the computed threshold value for each cell.
      feature.threshold <- rep(threshold.value, length(cell.id))
      # Create a boolean vector indicating whether each cell's
      # MALAT1 expression exceeds the computed threshold.
      feature.pass <- feature.logcounts[[feature]] > feature.threshold
      # Combine the cell IDs, threshold values, and pass/fail results into a data frame.
      results.df <- data.frame(
        cell.id = cell.id,
        feature.threshold = feature.threshold,
        feature.pass = feature.pass
      )

      # Rename the columns of the results data frame to include the feature name for clarity.
      names(results.df)[2:3] <- c(
        paste0(feature, ".threshold"),
        paste0(feature, ".pass")
      )

      return(results.df)
    }
  ) |> # Combine the results from all layers into a single data frame.
    data.table::rbindlist() |>
    as.data.frame()

  # Set the row names of the feature.metadata data frame to the cell IDs
  # for proper alignment with the Seurat object's metadata.
  row.names(feature.metadata) <- feature.metadata$cell.id
  # Delete the cell.id column from feature.metadata as it is now redundant with the row names.
  feature.metadata$cell.id <- NULL

  # Add the computed MALAT1 threshold and pass/fail results to the Seurat object's metadata.
  SeuratObject <- Seurat::AddMetaData(
    object = SeuratObject,
    metadata = feature.metadata
  )

  return(SeuratObject)
}
