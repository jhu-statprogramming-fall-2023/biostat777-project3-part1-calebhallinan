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
#' @export read_data_STTileMaker
#' @export STTileMaker
#' @export plot_tile
#'
#' @examples
#' # Make Tiles of each Spot in a Spatial Transcriptomics Dataset
#' directory = here("data/f12hr_140_processed/outs/")
#' tiles = STTileMaker_make_tiles(directory)


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


# read_data func
STTileMaker_make_tiles = function(data) {

  print("Making tiles...")
  # get positions of each spot
  pos.image = data$pos[rownames(data$posClean),]*data$scalefactor
  ## repeat for all tiles
  tiles = lapply(1:nrow(data$posClean), function(i) {
    # print(paste0("Making tile: ",i))
    ## for grabbing matrix elements in image
    pos.image.test = pos.image[i,]
    # get min and max of each tile
    pos.test.min = floor(pos.image.test - data$r)
    pos.test.max = ceiling(pos.image.test + data$r)

    # xs = seq(as.numeric(pos.test.min[1]), as.numeric(pos.test.max[1]))
    # ys = seq(as.numeric(pos.test.min[2]), as.numeric(pos.test.max[2]))
    xs = seq(as.numeric(pos.test.min[1]), as.numeric(pos.test.min[1])+23) # NOTE: need to fix this ?????
    ys = seq(as.numeric(pos.test.min[2]), as.numeric(pos.test.min[2])+23)

    # get image spot
    xspot = data$img[xs, ys, 1:3]

    return(xspot)
  })

  names(tiles) = rownames(data$posClean)
  print("Done!")

  return(tiles)

}

STTileMaker_plot_tile = function(tiles, index = 1) {
  # plot
  plot(grid::rasterGrob(tiles[[index]],
                        interpolate = FALSE,
                        width = grid::unit(1, "npc"),
                        height = grid::unit(1, "npc"))$raster)
}





# tiles = STTileMaker(directory)
#
#
library(here)
library(Matrix)
library(png)
library(rhdf5)
library(ggplot2)
library(rjson)

directory = here("data/f12hr_140_processed/outs/")
tiles = STTileMaker_read_data(directory)

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
