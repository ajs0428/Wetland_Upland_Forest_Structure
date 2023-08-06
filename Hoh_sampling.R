library(terra)
library(sf)
library(spatialEco)
library(lidR)
setwd("/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/")

hoh_WIP <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/HOH/Hoh_WIP_Mask0_10_2022.tif")

hoh_WIPw <- hoh_WIP >= 0.5
hoh_WIPu <- hoh_WIP < 0.5

GEE_HOHTCCproj <- rast("AGB/HOH/GEE_hohTCC.tif_proj.tif")

hoh_WIPwf <- mask(hoh_WIPw, GEE_HOHTCCproj >=50, maskvalues = 0, updatevalue = NA)
m <- c(FALSE, 0,
       TRUE, 1)
rclmat <- matrix(m, ncol=2, byrow=TRUE)
rc1 <- classify(hoh_WIPwf, rclmat, include.lowest=TRUE)

agFunc <- function(x) {sum(x*c(rep(0,2380),1,rep(0,2380)))}

hoh_WIPwf_agg <- terra::aggregate(rc1, 69, fun = agFunc)
plot(hoh_WIPwf_agg)

# create sample points from severity
pnts <- st_as_sf(as.points(hoh_WIPwf_agg))
plot(pnts)

wet_pnts <- pnts[pnts$WET > 0.5, ]

upl_pnts <- pnts[pnts$WET < 0.5, ] |>
    subsample.distance( 257, d = 500)


all_pnts <- rbind(wet_pnts, upl_pnts)

all_pnts_buff <- sf::st_buffer(all_pnts, 20)

# wf_pts <- spatSample(hoh_WIPwf, 1000, method="random", replace=FALSE, na.rm=TRUE, 
#            as.raster=FALSE, as.df=FALSE, as.points=TRUE, values=TRUE, cells=FALSE, 
#            xy=TRUE, ext=NULL, warn=TRUE, weights=NULL,  exp=5, exhaustive=FALSE)


ctg <- readLAScatalog("Lidar/hoh_2012/laz/")
ctg_poly <- st_as_sf(ctg)
st_crs(ctg_poly) <- 2927
st_crs(ctg) <- 2927
all_pnts_prj <- st_transform(all_pnts, "EPSG:2927")

small_pnts_prj <- st_crop(all_pnts_prj, ctg_poly)

x <- sf::st_coordinates(small_pnts_prj)[,1]
y <- sf::st_coordinates(small_pnts_prj)[,2]

opt_output_files(ctg) <- "clip_{}"
loc_ctg = clip_roi(ctg, small_pnts_prj, radius = 10)

#circle <- clip_circle(ctg, x, y, 20)
