## example of reading in Visium image, aligning to spots, making tiles

# starting from init.R ----------------------------------------------------
library(here)
## my directory
# dir <- '~/OneDrive - Johns Hopkins/Data_Private/Luigi_PRN_Visium/PRN_vs_wt/PRN1_wt/outs/'
dir = here("data/f12hr_140_processed/outs/")


# ## read in low resolution image of tissue for now
# library(png)
# img <- readPNG(paste0(dir, 'spatial/tissue_hires_image.png'))
# ## now we can plot with the image
# x <- grid::rasterGrob(img,
#                       interpolate = FALSE,
#                       width = grid::unit(1, "npc"),
#                       height = grid::unit(1, "npc"))
# library(ggplot2)
# plt <- ggplot() +
#   annotation_custom(
#     grob = x,
#     xmin = 0,
#     xmax = ncol(x$raster),
#     ymin = 0,
#     ymax = nrow(x$raster)) +
#   coord_fixed(
#     xlim = c(0, ncol(x$raster)),
#     ylim = c(0, nrow(x$raster))) +
#   theme_void()
# plt
#
# # read in gene expression
# #barcodes	string	Barcode sequences and their corresponding GEM wells (e.g. AAACGGGCAGCTCGAC-1)
# #data	uint32	Nonzero UMI counts in column-major order
# #indices	uint32	Zero-based row index of corresponding element in data
# #indptr	uint32	Zero-based index into data / indices of the start of each column, i.e., the #data corresponding to each barcode sequence
# #shape	uint64	Tuple of (# rows, # columns) indicating the matrix dimensions
# library(rhdf5)
# file <- paste0(dir, 'filtered_feature_bc_matrix.h5')
# h5f <- H5Fopen(file)
# geneInfo <- h5f$matrix
# genes <- geneInfo$features$name
# spotIDs <- geneInfo$barcodes
# library(Matrix)
# geneCounts <- sparseMatrix(
#   dims = geneInfo$shape,
#   i = as.numeric(geneInfo$indices),
#   p = as.numeric(geneInfo$indptr),
#   x = as.numeric(geneInfo$data),
#   index1 = FALSE
# )
# rownames(geneCounts) <- genes
# colnames(geneCounts) <- spotIDs
# dim(geneCounts)
#
# # get all spot barcode IDs and their positional coordinates
# spotPos <- read.csv(paste0(dir, "spatial/tissue_positions_list.csv"), header=TRUE, row.names=1)
# head(spotPos)
# pos <- spotPos[spotPos$X0==1,4:5] # Note: only if in_tissue a column
# colnames(pos) <- c('x', 'y')
# rownames(pos) <- rownames(spotPos)[spotPos$X0==1]
# head(pos)
# plot(pos)
# geneCounts <- geneCounts[, rownames(pos)]
#
# ## from scalefactor_hiress_json.json for lowres image
# ## alignment https://github.com/JEFworks-Lab/STdeconvolve/issues/39
# #scalefactor_hires <- 0.05881776 #lowres
# library(rjson)
# scale_factors = rjson::fromJSON(file=paste0(dir, 'spatial/scalefactors_json.json'))
# scalefactor_hires<- scale_factors$tissue_hires_scalef
#
# posAligned <- data.frame(x = pos[,2], y = -pos[,1] + dim(img)[1]/scalefactor_hires)*scalefactor_hires
# rownames(posAligned) <- rownames(pos)
# plot(posAligned,
#      xlim=c(-dim(img)[1], dim(img)[1]),
#      ylim=c(-dim(img)[2], dim(img)[2]))
# abline(h = 0)
# abline(v = 0)
#
# ## spot size
# r = 55*scalefactor_hires ## microns to pixels
#
# ## plot a gene
# ## gene and spot pos must be in same order
# head(sort(rowSums(geneCounts), decreasing=TRUE))
# g <- 'Spp1'
# ggexp <- scale(log10(geneCounts[g,rownames(posAligned)]+1))[,1]
# df <- data.frame(posAligned, ggexp)
# plt + geom_point(aes(x = x, y = y,
#                      color = ggexp),
#                  size = r/5,
#                  data = df) + ggtitle(g) +
#   scale_colour_gradient2(
#     low = "blue",
#     mid = "white",
#     high = "red")
#
#
# # clean data --------------------------------------------------------------
#
# ## arbitrary thresholds for now
# library(Matrix)
# table(colSums(geneCounts) > 1000)
#
# vi <- colnames(geneCounts)[colSums(geneCounts) > 1000]
# posClean <- posAligned[vi,]
# geneCountsClean <- geneCounts[, rownames(posClean)]
#
# g <- 'Spp1'
# ggexp <- scale(log10(geneCountsClean[g,rownames(posClean)]+1))[,1]
# df <- data.frame(posClean, ggexp)
# plt + geom_point(aes(x = x, y = y,
#                      color = ggexp),
#                  size = r/5,
#                  data = df) + ggtitle(g) +
#   scale_colour_gradient2(
#     low = "blue",
#     mid = "white",
#     high = "red")
#
# # make tiles --------------------------------------------------------------
# pos.image <- pos[rownames(posClean),]*scalefactor_hires
#
# ## given a spot position, get corresponding image region
# i <- 100
# pos.test <- posClean[i,] ## for plotting
# pos.image.test <- pos.image[i,] ## for grabbing matrix elements in image
# pos.test
# pos.image.test
# dim(img)
#
# pos.test.min <- floor(pos.image.test - r)
# pos.test.max <- ceiling(pos.image.test + r)
#
# xs <- seq(as.numeric(pos.test.min[1]), as.numeric(pos.test.max[1]))
# ys <- seq(as.numeric(pos.test.min[2]), as.numeric(pos.test.max[2]))
#
# xspot <- grid::rasterGrob(img[xs, ys, 1:3],
#                           interpolate = FALSE,
#                           width = grid::unit(1, "npc"),
#                           height = grid::unit(1, "npc"))
# plot(xspot$raster)
#
# pltspot <- ggplot() + geom_point(data=as.data.frame(posClean), aes(x=x, y=y), col='lightgrey') +
#   annotation_custom(
#     grob = xspot,
#     xmin = as.numeric(pos.test[1]-r),
#     xmax = as.numeric(pos.test[1]+r),
#     ymin = as.numeric(pos.test[2]-r),
#     ymax = as.numeric(pos.test[2]+r)) +
#   coord_fixed(
#     xlim = c(0, ncol(img)),
#     ylim = c(0, nrow(img)))
# pltspot
#
# ## repeat for all tiles
# tiles <- lapply(1:nrow(posClean), function(i) {
#   print(i)
#   pos.image.test <- pos.image[i,] ## for grabbing matrix elements in image
#
#   pos.test.min <- floor(pos.image.test - r)
#   pos.test.max <- ceiling(pos.image.test + r)
#
#   #xs <- seq(as.numeric(pos.test.min[1]), as.numeric(pos.test.max[1]))
#   #ys <- seq(as.numeric(pos.test.min[2]), as.numeric(pos.test.max[2]))
#   xs <- seq(as.numeric(pos.test.min[1]), as.numeric(pos.test.min[1])+23)
#   ys <- seq(as.numeric(pos.test.min[2]), as.numeric(pos.test.min[2])+23)
#
#   xspot <- img[xs, ys, 1:3]
#
#   return(xspot)
# })
# names(tiles) <- rownames(posClean)
# save(tiles, file=here("data/f12hr_140_processed/f12hr_140_processed_tiles.RData"))
#
# ## double check we are grabbing the right tiles
# ## (not very efficient...runs out of memory for all tiles but you get the sense)
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
#
# ## note to self, there are two coordinate systems here
# ## one image coordinates for grabbing pixels
# ## another of cartesian coordinates for plotting
#
#
# # simple check of rgb intensities -----------------------------------------
#
# ## given all these tiles, let's compute a "feature"
# ## just rum of rgb intensities
# ## would expect "off tissue" to be approx white
# plot(grid::rasterGrob(tiles[[1]],
#                       interpolate = FALSE,
#                       width = grid::unit(1, "npc"),
#                       height = grid::unit(1, "npc"))$raster)
# plot(grid::rasterGrob(tiles[[10]],
#                       interpolate = FALSE,
#                       width = grid::unit(1, "npc"),
#                       height = grid::unit(1, "npc"))$raster)
#
# 1 - mean(tiles[[1]])
# 1 - mean(tiles[[10]])
#
# features <- sapply(1:length(tiles), function(i) {
#   1 - mean(tiles[[i]])
# })
# names(features) <- names(tiles)
# features
#
# ## check putative empties
# names(which(features < 0.31))
#
# pltspot <- ggplot() + geom_point(data=as.data.frame(posClean), aes(x=x, y=y), col='lightgrey')
# pltspot
# xspotplots <- lapply(names(which(features < 0.31)), function(xspotname) {
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
#
# ## hypothesize this to correlate with total gene expression
# totgexp <- colSums(geneCountsClean)
# hist(log10(totgexp+1))
# hist(log10(totgexp[names(which(features < 0.31))]+1), add=TRUE) ## can still be oddly high? diffusion?
#
# plot(features, log10(totgexp+1)) ## no obvious correlation...
#
# ## what about one particular gene
# totgexp <- geneCountsClean[g,]
# hist(log10(totgexp+1))
# hist(log10(totgexp[names(which(features < 0.31))]+1), add=TRUE) ## can still be oddly high? diffusion?
#
# plot(features, log10(totgexp+1)) ## no obvious correlation...
#
# ## what about the other way around, look for spots with fewer genes and see what they look like
# totgexp <- geneCountsClean[g,]
# range(totgexp)
# table(totgexp < 10)
# names(which(totgexp < 10))
#
# pltspot <- ggplot() + geom_point(data=as.data.frame(posClean), aes(x=x, y=y), col='lightgrey')
# pltspot
# xspotplots <- lapply(names(which(totgexp < 10)), function(xspotname) {
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
#
# ## definitely grabbing right spots
#
#
# # do any genes correlate with simply rgb features -------------------------
#
# ## features need to have some variation
# features <- sapply(1:length(tiles), function(i) {
#   #1 - mean(tiles[[i]])
#   #mean(tiles[[i]][,,1]) ## red
#   mean(tiles[[i]][,,1] - tiles[[i]][,,3]) ## red minus blue
# })
# names(features) <- names(tiles)
# hist(features)
#
# ## normalize
# normFactor <- colSums(geneCountsClean)
# mat <- Matrix::t(Matrix::t(geneCountsClean)/normFactor)
# mat <- log10(mat * 1e6 +1)
#
# ## simple correlation
# plot(mat[g,], features)
# cor(mat[g,], features)
#
# cors <- sapply(rownames(mat), function(g) {
#   cor(mat[g,], features)
# })
# hist(cors)
# range(cors, na.rm=TRUE)
#
# names(which(cors == min(cors, na.rm=TRUE)))
# names(which(cors == max(cors, na.rm=TRUE)))
#
# g <- "Fth1"
# g <- 'Ehf'
# g <- "Retnla"
# ggexp <- scale(mat[g,])[,1]
# df <- data.frame(posClean, ggexp)
# plt + geom_point(aes(x = x, y = y,
#                      color = ggexp),
#                  size = r/5,
#                  data = df) + ggtitle(g) +
#   scale_colour_gradient2(
#     low = "blue",
#     mid = "white",
#     high = "red")
#
# plot(mat[g,], features)
# cor(mat[g,], features)
#
# ## note since I'm using the low res image, tiles are very small, not many features beyond intensities
# dim(tiles[[1]])
