library(terra)
library(sf)
library(spatialEco)
library(lidR)
library(tidyterra)
library(dplyr)
setwd("/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/")

hoh_WIP <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/HOH/Hoh_WIP_Mask0_10_2022.tif")

hoh_WIPw <- hoh_WIP >= 0.5
hoh_WIPu <- hoh_WIP < 0.5

hoh_poly <- vect("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/HOH/HOH_POLYGON_7_11_2022/")
GEE_HOHTCCproj <- rast("AGB/HOH/GEE_hohTCC.tif_proj.tif")
rip <- vect("NHD_surfacewater/NHD_H_17100101_HU8_Shape/Shape/NHDFlowline.shp") |> project(hoh_poly)
hoh_rip <- crop(rip, hoh_poly)
hoh_rip_fil <- vect("Wetland_Upland_Forest_Structure/data/hoh_rip_fil")#hoh_rip |> tidyterra::filter(visibility > 0) |> terra::buffer(15.24)

hoh_WIPwf <- hoh_WIPw |> mask(GEE_HOHTCCproj >=50, maskvalues = 0, updatevalue = NA) |>
    mask(hoh_rip_fil, inverse = T, updatevalue = NA)
hoh_WIPwf_2927 <- terra::project(hoh_WIPwf, "EPSG:2927")
plot(hoh_WIPwf, main = "hoh_WIPwf - wet forest")



m <- c(FALSE, 0,
       TRUE, 1)
rclmat <- matrix(m, ncol=2, byrow=TRUE)
rc1 <- classify(hoh_WIPwf, rclmat, include.lowest=TRUE)
plot(rc1)

agFunc <- function(x) {sum(x*c(rep(0,2380),1,rep(0,2380)))}

hoh_WIPwf_agg <- terra::aggregate(rc1, 69, fun = agFunc, filename = "Wetland_Upland_Forest_Structure/data/hoh_wip_binary_agg.tif",
                                  overwrite = T) # Aggregate to coarser resolution to avoid spatial autocorrelation
plot(hoh_WIPwf_agg)

# create sample points from aggregated raster
pnts <- st_as_sf(as.points(hoh_WIPwf_agg))
plot(pnts)

wet_pnts <- pnts[pnts$WET > 0.5, ]

upl_pnts <- pnts[pnts$WET < 0.5, ] |>
    subsample.distance( 108, d = 500)


all_pnts <- rbind(wet_pnts, upl_pnts)
plot(all_pnts)

wet_pnts_prj <- st_transform(wet_pnts, "EPSG:2927")
upl_pnts_prj <- st_transform(upl_pnts, "EPSG:2927")
all_pnts_prj <- st_transform(all_pnts, "EPSG:2927")
all_pnts_prj <- all_pnts_prj |> mutate(wetup = case_when(WET == 1 ~ "WET",
                                                         WET == 0 ~ "UPL",
                                                         .default = "other"))

plot(all_pnts_prj["wetup"])
all_pnts_prj_4326 <- st_transform(all_pnts_prj, "EPSG:4326")
writeVector(vect(all_pnts_prj_4326), "Wetland_Upland_Forest_Structure/data/all_pnts_prj_4326/all_pnts_prj_4326.shp",
            overwrite = T)

all_pnts_buff <- all_pnts |> sf::st_buffer(20)
plot(all_pnts_buff)

wf_pts <- spatSample(hoh_WIPwf, 1000, method="random", replace=FALSE, na.rm=TRUE,
           as.raster=FALSE, as.df=FALSE, as.points=TRUE, values=TRUE, cells=FALSE,
           xy=TRUE, ext=NULL, warn=TRUE, weights=NULL,  exp=5, exhaustive=FALSE)
plot(wf_pts)

ctg <- readLAScatalog("Lidar/hoh_all_pcs/datasetsA/")
plot(ctg)
ctg_poly <- ctg |> st_as_sf() |> sf::st_union()
st_crs(ctg_poly) <- 2927
st_crs(ctg) <- 2927

opt_output_files(ctg) <- "/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/Lidar/hoh_all_pcs/processing/hoh_{XCENTER}_{YCENTER}"
opt_laz_compression(ctg) <- TRUE
opt_merge(ctg) = FALSE

all_ctg <- clip_roi(ctg, all_pnts_prj, radius = 49.2) #15m radius, 30m pixel
#wet_ctg <- clip_roi(ctg, wet_pnts_prj, radius = 49.2) #15m radius, 30m pixel
#upl_ctg <- clip_roi(ctg, upl_pnts_prj, radius = 49.2) #15m radius, 30m pixel
(readLAS(all_ctg$filename[70])) #test

pc_filelist <- list.files("Lidar/hoh_all_pcs/processing/", pattern = ".laz")

pc_pt_func <- function(filelist){
    pc_mat <- matrix(nrow = length(filelist), ncol = 3)
    for (i in 1:length(filelist)){
        pc_mat[i, 1] = i
        pc_mat[i, 2] = stringr::str_extract(filelist[i], "(?<=_)[^_]+(?=_)")
        pc_mat[i, 3] = stringr::str_extract(filelist[i], "(?<=_)[^_]+(?=.laz)")
    }
    pc_dat <- as.data.frame(pc_mat)
    colnames(pc_dat) <- c("n", "x", "y")
    pc_dat <- pc_dat |> mutate(x = as.numeric(x),
                                y = as.numeric(y))
    pc_vect <- vect(pc_dat, geom = c("x", "y"), crs = "EPSG:2927", keepgeom = TRUE)
    return(pc_vect)
}

hoh_pc_pts <- pc_pt_func(pc_filelist) |> 
    terra::extract(x = hoh_WIPwf_2927, bind = TRUE) |> 
    mutate(wetup = case_when(WET >= 0.5 ~ "WET",
                             WET < 0.5 ~ "UPL",
                             .default = "other"))

wetupl_table <- data.frame(pc_filelist) |> cbind(hoh_pc_pts$wetup) |> 
    dplyr::rename(wetup = `hoh_pc_pts$wetup`)
wetupl_table |> group_by(wetup) |> tally()

#---------------------------------------------------------------------------------------
wet_list <- list()
for (i in c(1:length(wet_ctg))){
    #print(i)
    wet_list[[i]] <- normalize_height(wet_ctg[[i]], knnidw())
}
plot(wet_list[[1]], color = "Classification")
upl_list <- list()
for (i in c(1:length(upl_ctg))){
    #print(i)
    upl_list[[i]] <- normalize_height(upl_ctg[[i]], knnidw())
}

wet_metrics <- cloud_metrics(wet_ctg[[1]], func = .stdmetrics_z)
wet_metrics_df <- as.data.frame(wet_metrics)
upl_metrics <- cloud_metrics(upl_ctg[[1]], func = .stdmetrics_z)
upl_metrics_df <- as.data.frame(upl_metrics)
# wet_metrics2 <- cloud_metrics(wet_ctg[[2]], func = .stdmetrics_z)
# wet_metrics2_df <- as.data.frame(wet_metrics2)
# test <- rbind(wet_metrics_df, wet_metrics2_df)

columns <- names(wet_metrics_df)
wetData = data.frame(matrix(nrow = 0, ncol = length(columns)))
names(wetData) <- names(wet_metrics_df)
uplData = data.frame(matrix(nrow = 0, ncol = length(columns)))
names(uplData) <- names(wet_metrics_df)

for(i in 1:length(wet_list)){
    wetData[i,] <- cloud_metrics(wet_list[[i]], func = .stdmetrics_z)
    #wet_mets_df <- as.data.frame(wet_mets)
}
wetData$wetupl <- "WET"

for(i in 1:length(upl_list)){
    uplData[i,] <- cloud_metrics(upl_list[[i]], func = .stdmetrics_z)
    #upl_mets_df <- as.data.frame(upl_mets)
}
uplData$wetupl <- "UPL"

myData <- rbind(wetData, uplData)
myData_num <- (myData[, c(-6,-37)])


myData_clean <- 
pc <- prcomp((myData[, c(-6,-37)]),  
                         scale = TRUE)
PCA1 <- pc$x[,1]
PCA2 <- pc$x[,2]


library(ggplot2)
ggplot((myData[, c(-6,-37)]), 
       aes(x = PCA1, 
           y =  PCA2, 
           color = myData[, 37])) +
    geom_point() +
    stat_ellipse()


library(ggfortify)
autoplot(pc, data = (myData[, c(-6)]), colour = "wetupl", loadings = TRUE, loadings.label = TRUE)


upl_metrics <- cloud_metrics(upl_ctg, func = .stdmetrics_z) 



#######

test <- readLAS("Lidar/hoh_all_pcs/processing/hoh_799099.3_940639.6.laz")
plot(test)
