#' Plot a tile of a spot in a Spatial Transcriptomics Dataset
#'
#' Plot a tile of a spot in a Spatial Transcriptomics Dataset
#'
#' @details This function will make plot a tile of a spot in a Spatial Transcriptomics Dataset as a quality control check
#'
#' @param tiles, tiles object from function STTileMaker_make_tiles
#' @param index, the index of a tile you want to visualize, range is 1 to n number of tiles
#'
#' @return Plot of a single tile
#'
#' @import here
#' @import Matrix
#' @importFrom png readPNG
#' @import rhdf5
#' @import ggplot2
#' @import rjson
#' @export STTileMaker_plot_tile
#'
#' @examples
#' STTileMaker_plot_tile(tiles, 10)


STTileMaker_plot_tile = function(tiles, index = 1) {
  # plot
  plot(grid::rasterGrob(tiles[[index]],
                        interpolate = FALSE,
                        width = grid::unit(1, "npc"),
                        height = grid::unit(1, "npc"))$raster)
}


