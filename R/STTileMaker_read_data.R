#' Reads in the data that contains the image, gene expression, and various other variables
#'
#' Reads in the data that contains the image, gene expression, and various other variables
#'
#' @details This function will reads in the data that contains the image, gene expression, and various other variables
#'
#' @param directory, directory where the ST data is located
#'
#' @return list of lists with:
#'
#' pos - original positions of spots
#'
#' posAligned - positions of spots after data alignment
#'
#' posClean - positions of spots after data cleaning
#'
#' img - the H&E image in matrix form
#'
#' geneCounts - matrix of the of the raw data
#'
#' geneCountsClean - matrix that has been cleaned with standard ST data methods
#'
#' r - scale for resizing spots
#'
#' scalefactor - scalefactor from image to cartesian scales
#'
#' @import here
#' @import Matrix
#' @importFrom png readPNG
#' @import rhdf5
#' @import ggplot2
#' @import rjson
#' @export STTileMaker_read_data
#'
#' @examples
#' data = STTileMaker_read_data(here::here("data/f12hr_140_processed/outs/"))



# read_data func
STTileMaker_read_data <- function(directory) {
  # stop if it is not a character object
  stopifnot(is.character(directory))
  # read in image
  img = readPNG(paste0(directory, 'spatial/tissue_hires_image.png'))

  # file name
  gene_data_file = paste0(directory, 'filtered_feature_bc_matrix.h5')
  # read in file
  # h5f = H5Fopen(gene_data_file) # NOTE: this method was not letting me knit documents, saying file was unreadable
  # get gene expression (list of objects)
  # geneInfo = h5f$matrix

  # get gene expression (list of objects)
  geneInfo <- h5read(gene_data_file, "/matrix")
  # get gene names (vector of genes)
  genes = geneInfo$features$name
  # get spot ids (vector of names)
  spotIDs = geneInfo$barcodes
  # make sparse matrix of expression
  geneCounts = sparseMatrix(
    # set dimensions
    dims = geneInfo$shape,
    # set indices
    i = as.numeric(geneInfo$indices),
    # make numeric
    p = as.numeric(geneInfo$indptr),
    # add gene expression
    x = as.numeric(geneInfo$data),
    index1 = FALSE
  )
  # set rownames to gene names
  rownames(geneCounts) = genes
  # set col names to spot ids
  colnames(geneCounts) = spotIDs
  # print dim
  print(paste0("Gene Matrix Dimensions: ", dim(geneCounts)))

  # read in spot positions
  spotPos = read.csv(paste0(directory, "spatial/tissue_positions_list.csv"), header=TRUE, row.names=1)
  # grab x and y of spots
  pos = spotPos[spotPos$X0==1,4:5]
  # rename to x and y
  colnames(pos) = c('x', 'y')
  # grab only spots that are in the tissue
  rownames(pos) = rownames(spotPos)[spotPos$X0==1] # Note: hopefully always X0?
  print(head(pos))
  # get genecounts of only spots in tissue
  geneCounts = geneCounts[, rownames(pos)]

  # read in scalefactor
  scalefactors_rjson = rjson::fromJSON(file=paste0(directory, 'spatial/scalefactors_json.json'))
  # get hires
  scalefactor = scalefactors_rjson$tissue_hires_scalef
  # align positions
  posAligned = data.frame(x = pos[,2], y = -pos[,1] + dim(img)[1]/scalefactor)*scalefactor
  rownames(posAligned) = rownames(pos)

  # spot size
  r = 55*scalefactor ## microns to pixels

  ## clean data

  ## arbitrary thresholds for now
  table(colSums(geneCounts) > 1000)
  # get count above threshold
  vi = colnames(geneCounts)[colSums(geneCounts) > 1000]
  posClean = posAligned[vi,]
  geneCountsClean = geneCounts[, rownames(posClean)]

  return_object = list("geneCountsClean" = geneCountsClean, "posClean" = posClean)

  # make everything to list
  return_object = list("img" = img, "geneCounts" = geneCounts, "pos" = pos, "posAligned" = posAligned, "geneCountsClean" = geneCountsClean, "posClean" = posClean, "r" = r, "scalefactor" = scalefactor)

  # return
  return (return_object)
}



