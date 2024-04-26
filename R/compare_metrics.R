library(dplyr)
library(readr)
library(ggplot2)
library(foreach)
library(lidR)

source("R/metrics.R")

clouds <- read_csv("data/wetupl_table_for_ptclouds.csv")

metric_results <- lapply(
  clouds$pc_filelist, 
  function(f) {
    this_cloud <- readLAS(file.path("data", "processing", f))
    
    c(
      # Remember what the file is called
      list(pc_filelist=f),
      # Standard metrics in lidR
      cloud_metrics(this_cloud, ~stdmetrics_z(Z)),
      # Other stuff
      list(
        crr=canopy_relief_ratio(this_cloud),
        fhd=foliar_height_diversity(this_cloud),
        fd=fractal_dimension(this_cloud),
        fds=fractal_dimension_shannon(this_cloud),
        och=outer_canopy_height(this_cloud),
        rug=rugosity(this_cloud),
        rum=rumple_index_las(this_cloud)
      )
    )
  }
)

metric_df <- bind_rows(metric_results) %>%
  right_join(clouds, by="pc_filelist")

# Apply scaling
metrics_scale <- metric_df %>%
  select(-wetup, -pc_filelist, -pzabove2) %>%
  as.matrix() %>%
  scale()

# 100% of first returns are above 2f in all clouds so ignore that column
metric_pca <- princomp(metrics_scale)

metric_df$PCA1 <- metric_pca$scores[, 1]
metric_df$PCA2 <- metric_pca$scores[, 2]

metric_df %>%
  filter(wetup != "other") %>%
  ggplot(aes(x=PCA1, y=PCA2, color=wetup)) +
  geom_point() + stat_ellipse() +
  theme_bw()

write_csv(metric_df, "data/metrics.csv")
