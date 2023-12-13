#' Make tiles of each spot in a Spatial Transcriptomics Dataset
#'
#' Make tiles of each spot in a Spatial Transcriptomics Dataset
#'
#' @details This function will make Tiles of each Spot in a Spatial Transcriptomics Dataset when given
#' the correct input
#'
#' @param data, object from function STTileMaker_read_data
#'
#' @return H&E tiles from each spot is returned as a list of lists
#'
#' @import here
#' @import Matrix
#' @importFrom png readPNG
#' @import rhdf5
#' @import ggplot2
#' @import rjson
#' @export STTileMaker_make_tiles
#'
#' @examples
#' tiles = STTileMaker_make_tiles(data)


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
