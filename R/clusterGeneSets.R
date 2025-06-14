#' Cluster gene sets
#'
#' @param geneSetsList 'gsList' object created with the 'initializeList()' function.
#' @param resolution Vector with the desired clustering resolutions to use. Higher number usually results in more clusters.
#' @param verbose a logic indicating whether to print progress statements. Defaults to TRUE.
#'
#' @return 'gsList' object with clustering stored. The best resolution (highest
#'  silhouette score) is annoated and used throughout the rest of the workflow (unless
#'  specified in subsequent functions).
#' @export
#'
#' @examples
#' # Create gsListObject
#' geneSetsList <- initializeList(geneSetsList = exampleGeneSets, geneSetsDF = exampleEnrichment)
#' gsListObject <- clusterGeneSets(geneSetsList = geneSetsList)
clusterGeneSets <-
  function(
    geneSetsList,
    resolution = c(0.1, 0.3, 0.6, 1.1, 1.6, 2.1, 3.1, 4.1, 5.1),
    verbose = TRUE
) {

    # Check input
    if (!class(geneSetsList) == "gsList") {
      stop("'geneSetsList' should be a gsList object.")
    }

    if (!is.numeric(resolution)) {
      stop("'resolution' should be numeric.")
    }

    if(verbose) {
      message("Clustering gene sets...")
    }

    # Create seurat object
    suppressWarnings(
      seuratObj <-
        Seurat::CreateSeuratObject(counts = geneSetsList@genesetsDataframe,
                                   assay = "GSEA")
    )

    seuratObj <- Signac::RunTFIDF(seuratObj,
                                  verbose = F,
                                  method = 3) # verified by new grid search

    seuratObj <- Signac::FindTopFeatures(seuratObj,
                                         min.cutoff = "q0")

    seuratObj <- suppressWarnings(Signac::RunSVD(seuratObj,
                                                 approx = F,
                                                 verbose = F))

    # get max dim
    elbowRes <- kvsElbow(variance = seuratObj@reductions$lsi@stdev)

    seuratObj <- suppressWarnings(
      Seurat::RunUMAP(
        seuratObj,
        reduction = "lsi",
        dims = 2:elbowRes$elbow_cutoff, # verified by new grid search
        verbose = F
      )
    )

    seuratObj <- Seurat::FindNeighbors(
      object = seuratObj,
      reduction = 'lsi',
      dims = 2:elbowRes$elbow_cutoff, # verified by new grid search
      annoy.metric = "cosine",        # verified by new grid search
      k.param = 30,                   # updated by new grid search
      compute.SNN = "True",
      verbose = F
    )

    ### Find clusters
    seuratObj2 <- NULL
    i <- length(resolution)
    while(is.null(seuratObj2)) {
      suppressWarnings(
        seuratObj2 <- tryCatch({
          Seurat::FindClusters(
            object = seuratObj,
            resolution = resolution[1:i],
            verbose = F
          )
        }, error = {
          function(x) {
            NULL
          }
        })
      )

      i <- i -1
      if(i < 1) {
        message("No valid resolution. Please try with other resolutions.")
        return(NULL)
      }

    }
    seuratObj <- seuratObj2


    ### Set optimal clustering as default
    if(TRUE) {
      # Extract silluette scores
      meanSil <- sapply(resolution[1:(i+1)], function(localRes) {

        localClust <- as.integer(as.character(seuratObj@meta.data[,stringr::str_c('GSEA_snn_res.',localRes)]))

        if(length(unique(localClust)) > 1) {
          silInfo <- base::summary(
            cluster::silhouette(
              x = localClust,
              dmatrix = 1 - seuratObj@graphs$GSEA_snn
            )
          )
          mean_clust_sil <- mean(silInfo$clus.avg.widths)
          return(mean_clust_sil)
        } else {
          return(0)
        }

      })

      ### Extract best resolution
      bestRes <- resolution[which(meanSil == max(meanSil))[1] ]

      ### Set as default
      Seurat::Idents(seuratObj) <- seuratObj@meta.data[,stringr::str_c('GSEA_snn_res.',bestRes)]
      seuratObj@meta.data$best_cluster_res <- stringr::str_c('GSEA_snn_res.',bestRes)

    }

    # Save seurat object to list
    geneSetsList@cluster <- seuratObj

    # Add meta data column with gene set names
    geneSetsList@cluster <-
      Seurat::AddMetaData(geneSetsList@cluster,
                          geneSetsList@cluster@assays[["GSEA"]]$counts@Dimnames[[2]],
                          col.name = "gene_set")

    if(verbose) {
      for (i in resolution) {
        resolutionColumn <- paste0("GSEA_snn_res.", i)
        clusterNumber <-
          length(levels(geneSetsList@cluster@meta.data[[resolutionColumn]]))

        message(paste0(clusterNumber, " clusters found for resolution ", i))
      }
    }

    if(verbose) {
      message(paste(
        'Determined',
        bestRes,
        "as the optimal clustering resolution. This is set as default.",
        sep = ' '
      ))
    }

    if(verbose) {
      message("Clustering DONE")
    }

    return(geneSetsList)
  }
