#### Landsat 1-5 NDVI Processing & Analysis ####
# Kelsey E. Carter
# March 2018 (Revised September 2019)

#### Script Purpose ####
# The purpose of the script is to better understand spatial patterns of plant growth over time. The script has two parts: Part I processes Landsat imagery for NDVI and Part II documents exploratory analyses and visualizations of NDVI over time. More details are provided in each section.

##### load libraries #####
library(raster)
library(rgdal)
library(rgeos)
library(dplyr)
library(maptools)
library(png)
library(sp)
library(knitr)
library(RColorBrewer)
library(rasterVis)
library(ggplot2)
library(jpeg)
library(lubridate)
# library(googledrive)

#### Helpful functions ####

# `norm_diff`:  The `norm_diff` function finds the difference between two objects over the sum of the two objects (normalized difference).
# Input: band 1, band 2, where band 1 will be subtracted from band 2 (band2 - band1). The input will be a raster. For NDVI, b1 = red band and b2 = NIR band
# Output: A raster file of the normalized difference
norm_diff <- function(b1, b2) {
  (b2 - b1)/(b1 + b2)
}

# The `check_create_dir` function checks to see if a folder exists in the directory. 
# Input: directory path that you want to check/create
# Output: the creation of the directory and/or a statement of the action taken by the function
check_create_dir <- function(dir_path) {
  if(dir.exists(dir_path)) { # if directory DOES exist...
    print("the directory already exists") # ... print this statement
  }
  if(!dir.exists(dir_path)) { # if directory does NOT exist (!)...
    dir.create(dir_path) # ... create the directory and...
    print("the directory has been created") # ...print this statement
  }
}

#### Part I: Landsat processing of NDVI ####
# This part of the script processes Landsat satellite imagery from Landsat 1-5 (1972-2013) to calculate the normalized difference vegetation index (NDVI), which is a metric of plant growth and health. Images were downloaded prior to running the script and are publicly available via the USGS EarthExplorer at https://earthexplorer.usgs.gov/ (Landsat 1-5 MSS Collection 1 Level 1). The MSS sensor was used because  the research question required multiple decades of data and MSS is the longest running sensor. Only images that were 100% cloud-free within the study site boundary were included.

##### Load in the spatial extent file #####
sensor_boundary <- readOGR("data/study_area_boundary/sensor_polygon_complex2.shp")

# Get list of landsat 1-5 zipped files:
landsat_zip_filelist <- list.files("data/landsat_1-5_MSS/",
                                   full.names = TRUE)

# Create and write blank dataframe
landsat_df <- data.frame(landsat_mission = as.numeric(),
                         sensor = as.character(),
                         date = as.numeric(),
                         pixel_count = as.numeric(),
                         ndvi_max = as.numeric(),
                         ndvi_min = as.numeric(),
                         ndvi_mean = as.numeric(),
                         ndvi_sd = as.numeric(),
                         ndvi_cv = as.numeric())
output_csv_name <- "output/landsat_ndvi.csv"
check_create_dir("output/")
write.csv(x = landsat_df, 
          file = output_csv_name,
          row.names = FALSE)

# Create folder for unzipped Landsat data files
check_create_dir("data/landsat_1-5_MSS_unzipped")

##### For loop of landsat data #####
for (i in seq_along(landsat_zip_filelist)) {
  
  # Define location of files
  unzipped_directory <- paste0("data/landsat_1-5_MSS_unzipped/", basename(landsat_zip_filelist[i]))
  
  # Unzip the file from the tarball (.tar.gz)
  unzip_file <- untar(landsat_zip_filelist[i], 
                      list = FALSE,
                      exdir = unzipped_directory)
  
  # Get list of .TIF files in new unzipped folder
  unzip_bands_list <- list.files(path = unzipped_directory, 
                                 pattern = glob2rx(pattern = "*_B*.TIF$"),
                                 full.names = TRUE,
                                 recursive = FALSE)
  
  # Create stack, then brick of the bands
  landsat_st <- stack(unzip_bands_list)
  landsat_br <- brick(landsat_st)
  
  # convert crs of sensor boundary to the crs of the sensor_boundary
  sensor_boundary_utm <- spTransform(sensor_boundary, CRS = crs(landsat_br))
  
  # Crop landsat brick
  landsat_crop <- crop(landsat_br, sensor_boundary_utm)
  landsat_mask <- mask(landsat_crop, sensor_boundary_utm)
  
  # Calculate NDVI for landsat_mask (for Landsat 1-5)
  mss_ndvi <- overlay(landsat_mask[[2]], # red band (0.6-0.7 nm): band 5 (Landsat 1-3) or band 2 (Landsat 4-5)
                      landsat_mask[[4]], # NIR (0.8-1.1): band 7 (Landsat 1-3) or band 4(Landsat 4-5)
                      fun = norm_diff) 
  
  # Explore NDVI values for layer:
  scene_ndvi_max <- cellStats(mss_ndvi, stat = "max", na.rm = TRUE)
  scene_ndvi_min <- cellStats(mss_ndvi, stat = "min", na.rm = TRUE)
  scene_ndvi_mean <- cellStats(mss_ndvi, stat = "mean", na.rm = TRUE)
  scene_ndvi_sd <- cellStats(mss_ndvi, stat = "sd", na.rm = TRUE, asSample = TRUE)
  scene_ndvi_cv <- 100*(scene_ndvi_sd/scene_ndvi_mean)
  
  # Calculate number of pixels
  ndvi_pixels_extent <- mss_ndvi@data@values
  ndvi_pixels_crop <- na.omit(ndvi_pixels_extent) 
  ndvi_pixels_crop <- as.numeric(ndvi_pixels_crop)
  ndvi_pixel_count <- length(ndvi_pixels_crop)
  
  # Extract satellite/sensor information
  base_name_with_extension <- basename(landsat_zip_filelist[i])
  base_name <- substr(base_name_with_extension, 1, 21)
  landsat_sat <- substr(base_name, 3, 3)
  landsat_sensor <- substr(base_name, 2, 2)
  
  # Extract date
  scene_year <- as.numeric(substr(base_name, 10, 13))
  scene_doy <- as.numeric(substr(base_name, 14, 16))
  origin_date <- paste0((scene_year - 1), "-12-31")
  scene_date <- as.Date(scene_doy, origin = origin_date)
  scene_month <- month(scene_date)
  scene_day <- day(scene_date)
  date_as_character <- paste0()
  
  
  # Create dataframe of ndvi summary data
  temp_df <- data.frame(landsat_mission = landsat_sat,
                        sensor = landsat_sensor,
                        date = scene_date,
                        pixel_count = ndvi_pixel_count,
                        ndvi_max = scene_ndvi_max,
                        ndvi_min = scene_ndvi_min,
                        ndvi_mean = scene_ndvi_mean,
                        ndvi_sd = scene_ndvi_sd,
                        ndvi_cv = scene_ndvi_cv)
  
  # Load and merge with existing NDVI csv file
  old_csv <- read.csv("output/landsat_ndvi.csv") %>% 
    mutate(date = as.Date(date))
  merged_df <- rbind(old_csv, temp_df)
  write.csv(x = merged_df, 
            file = "output/landsat_ndvi.csv",
            row.names = FALSE)
  
  # Plot scene NDVI
  plot_ndvi <- plot(round(mss_ndvi, digits = 2),
                    main = scene_date,
                    # col = ndvi_colors,
                    col = rev(terrain.colors(11)),
                    breaks = seq(-1, 1, length.out = 11), 
                    pin = c(4,5),
                    legend = TRUE,
                    legend.args = list(text = 'NDVI Value', 
                                       side = 2, font = 2, line = 0, cex = 1))
  
  # Save NDVI plot to Landsat_NDVI_plots folder
  check_create_dir("output/NDVI_plots")
  dev.print(png, paste0("output/NDVI_plots/", "L", landsat_sat, "_", scene_date, ".png"), width = 500, height = 500)
  
  # Clear dev
  dev.new()
  
  # Print the progress...
  print(paste("Finished", i, " of ", length(landsat_zip_filelist)))
}
