#' Make Tiles of each Spot in a Spatial Transcriptomics Dataset
#'
#' Make Tiles of each Spot in a Spatial Transcriptomics Dataset
#'
#' @details This function will make Tiles of each Spot in a Spatial Transcriptomics Dataset when given
#' the correct input
#'
#' @param directory, directoryectory of ST data
#'
#' @return Tiles??
#'
#' @import here
#' @import Matrix
#' @importFrom png readPNG
#' @import rhdf5
#' @import ggplot2
#' @import rjson
#' @export STTileMaker
#'
#'
#' @examples
#' # Make Tiles of each Spot in a Spatial Transcriptomics Dataset
#' directory = here("data/f12hr_140_processed/outs/")
#' tiles = STTileMaker(directory)


# read_data func
read_data <- function(directory) {
  # stop if it is not an object
  stopifnot(is.character(directory))
    # read in image
  img = readPNG(paste0(directory, 'spatial/tissue_hires_image.png'))

  # file name
  file = paste0(directory, 'filtered_feature_bc_matrix.h5')
  # read in file
  h5f = H5Fopen(file)
  # get gene expression (list of objects)
  geneInfo = h5f$matrix
  # get gene names (vector of genes)
  genes = geneInfo$features$name
  # get spot ids (vector of names)
  spotIDs = geneInfo$barcodes
  # read in matrix
  library(Matrix)
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
  dim(geneCounts)


  # read in spot positions
  spotPos = read.csv(paste0(directory, "spatial/tissue_positions_list.csv"), header=TRUE, row.names=1)
  # head(spotPos)
  # grab x and y of spots
  pos = spotPos[spotPos$X0==1,4:5]
  # rename to x and y
  colnames(pos) = c('x', 'y')
  # grab only spots that are in the tissue
  rownames(pos) = rownames(spotPos)[spotPos$X0==1]
  print(head(pos))
  # plot
  # plot(pos)
  # get genecounts of only spots in tissue
  geneCounts = geneCounts[, rownames(pos)]

  # read in scalefactor
  scalefactors_rjson = rjson::fromJSON(file=paste0(directory, 'spatial/scalefactors_json.json'))
  # get hirest
  scalefactor = scalefactors_rjson$tissue_hires_scalef
  # align positions
  posAligned = data.frame(x = pos[,2], y = -pos[,1] + dim(img)[1]/scalefactor)*scalefactor
  rownames(posAligned) = rownames(pos)

  ## spot size
  r = 55*scalefactor ## microns to pixels = NEED FOR PLOTTING

  return_object = list("img" = img, "geneCounts" = geneCounts, "pos" = pos, "posAligned" = posAligned, "r" = r, "scalefactor" = scalefactor)

  return (return_object)
}



# clean_data func
clean_data = function(geneCounts, posAligned) {
  ## arbitrary thresholds for now
  table(colSums(geneCounts) > 1000)
  # get count above threshold
  vi = colnames(geneCounts)[colSums(geneCounts) > 1000]
  posClean = posAligned[vi,]
  geneCountsClean = geneCounts[, rownames(posClean)]

  return_object = list("geneCountsClean" = geneCountsClean, "posClean" = posClean)

  return(return_object)
}


# read_data func
STTileMaker = function(directory) {

  # read in data
  data = read_data(directory)
  # clean data
  data_clean = clean_data(data$geneCounts, data$posAligned)

  # get positions of each spot
  pos.image = data$pos[rownames(data_clean$posClean),]*data$scalefactor
  ## repeat for all tiles
  tiles = lapply(1:nrow(data_clean$posClean), function(i) {
    print(paste0("Making tile: ",i))
    ## for grabbing matrix elements in image
    pos.image.test = pos.image[i,]
    # get min and max of each tile
    pos.test.min = floor(pos.image.test - data$r)
    pos.test.max = ceiling(pos.image.test + data$r)

    #xs = seq(as.numeric(pos.test.min[1]), as.numeric(pos.test.max[1]))
    #ys = seq(as.numeric(pos.test.min[2]), as.numeric(pos.test.max[2]))
    xs = seq(as.numeric(pos.test.min[1]), as.numeric(pos.test.min[1])+23) # NOTE: need to fix this ?????
    ys = seq(as.numeric(pos.test.min[2]), as.numeric(pos.test.min[2])+23)

    # get image spot
    xspot = data$img[xs, ys, 1:3]

    return(xspot)
  })

  names(tiles) = rownames(data_clean$posClean)

  return(tiles)

}

# tiles = STTileMaker(directory)
#
#
# library(here)
# library(Matrix)
# library(png)
# library(rhdf5)
# library(ggplot2)
# library(rjson)
#
# directory = here("data/f12hr_140_processed/outs/")
# tiles = STTileMaker(directory)
#
# z = read_data(directory)
#
# z_clean = clean_data(z$geneCounts, z$posAligned)
#
#
# pltspot <- ggplot() + geom_point(data=as.data.frame(posClean), aes(x=x, y=y), col='lightgrey')
# pltspot
# xspotplots <- lapply(names(tiles)[1:100], function(xspotname) {
#
#   pos.test <- posClean[xspotname,]
#   pos.test.min <- floor(pos.test - r)
#   pos.test.max <- ceiling(pos.test + r)
#
#   xspot <- grid::rasterGrob(tiles[[xspotname]],
#                             interpolate = FALSE,
#                             width = grid::unit(1, "npc"),
#                             height = grid::unit(1, "npc"))
#
#   pltspot <<- pltspot + annotation_custom(
#     grob = xspot,
#     xmin = as.numeric(pos.test.min[1]),
#     xmax = as.numeric(pos.test.max[1]),
#     ymin = as.numeric(pos.test.min[2]),
#     ymax = as.numeric(pos.test.max[2]))
# })
# pltspot +
#   coord_fixed(
#     xlim = c(0, ncol(img)),
#     ylim = c(0, nrow(img)))
