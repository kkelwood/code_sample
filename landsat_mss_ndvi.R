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
library(ggplot2)
library(cowplot)
library(gridExtra)
library(plotrix)

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
# This part of the script processes Landsat satellite imagery from Landsat 1-5 (1972-2013) to calculate the normalized difference vegetation index (NDVI), which is a metric of plant growth and health. Images were downloaded prior to running the script and are publicly available via the USGS EarthExplorer at https://earthexplorer.usgs.gov/ (Landsat 1-5 MSS Collection 1 Level 1). The MSS sensor was used because  the research question required multiple decades of data and MSS is the longest running sensor. Only images that were 100% cloud-free within the study site boundary were included. For simplicity, only 5 satellite images are included in the accompanying data folder, but the original analysis used 54 images. Part II uses the analysis from all 54 images.

# Load in the spatial extent file
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

# For loop of landsat data
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
  
  # Save NDVI plot to Landsat_NDVI_maps folder
  check_create_dir("output/NDVI_maps")
  dev.print(png, paste0("output/NDVI_maps/", "L", landsat_sat, "_", scene_date, ".png"), width = 500, height = 500)
  
  # Clear dev
  dev.new()
  
  # Print the progress...
  print(paste("Finished", i, " of ", length(landsat_zip_filelist)))
}


#### Part II: Exploratory Data Analysis & Visualization ####

# Load in csv file with summary NDVI values (e.g. mean, max, min, sd, cv)
landsat_ndvi <- read.csv("output/landsat_ndvi_all.csv") %>% 
  mutate(date = as.POSIXct(date)) %>% 
  mutate(YEAR = year(date)) %>% 
  mutate(MONTH = month(date)) %>% 
  mutate(DAY = day(date)) %>% 
  mutate(MMDD = paste0(MONTH, "-", DAY)) %>% 
  mutate(MMDD = as.POSIXct(MMDD, format = "%m-%d")) %>% 
  mutate(DOY = yday(MMDD))

str(landsat_ndvi)

# How many years are represented?
length(unique(landsat_ndvi$YEAR))

# Which years are represented?
sort(unique(landsat_ndvi$YEAR))

# Histogram of how many scenes are available by week of the growing season
ggplot(landsat_ndvi, aes(x = week(MMDD))) +
  geom_histogram(binwidth = 0.5) + 
  labs(title = "How many cloud-free Landsat images are available?", subtitle = "Total from 1973-2013", x = "Week of the Year") + 
  theme_bw()

# Histogram of how many scenes are available each year
ggplot(landsat_ndvi, aes(x = YEAR)) +
  geom_histogram(binwidth = 0.5, 
                 center = 0) + 
  labs(title = "How many cloud-free images are available each year?",
       x = "Year") + 
  theme_bw()

# Histogram of NDVI max values
ggplot(landsat_ndvi, aes(x = ndvi_max)) +
  geom_histogram(binwidth = 0.1, 
                 center = 0) + 
    labs(title = "Histogram of max NDVI per scene",
         x = "NDVI Max") + 
  theme_bw()

# Create multi-panel plot of NDVI values (max, min, mean, sd, cv) over the growing season. Each point represents an MSS scene taken between 1973 and 1992.

## Define variables are of interest: mean, max, min, sd, cv
variable_list <- names(landsat_ndvi[5:9])

## Populate list `p` with a plot for each of the variables of interest
for (i in variable_list) {
  p[[i]] <- ggplot(data = landsat_ndvi, aes_string(x = "MMDD", y = i)) +
    geom_point(aes(color = YEAR)) +
    scale_color_gradientn(colors = terrain.colors(5), name = "Year") + 
    geom_smooth(color = "darkgrey") +
    labs(x = "Date", title = i) + 
    theme_bw()
  print(plot)
}

## Arrange all plots in list `p` in grid
do.call(grid.arrange, p)

# Create new dataframe that includes only overall max and min NDVI value for each year
landsat_ndvi_yearly_max <- landsat_ndvi %>% 
  group_by(YEAR) %>% 
  summarise(year_max = max(ndvi_max))
landsat_ndvi_yearly_max2 <- merge(landsat_ndvi, landsat_ndvi_yearly_max)

# Plot of max NDVI (of a single pixel) for each year
max_plot <- ggplot(data = landsat_ndvi_yearly_max2, aes(x = YEAR,  year_max)) +
  geom_point(aes(color = DOY)) +
  scale_color_gradient(low = "#a1dab4", high = "#253494", 
                       name = "Day of Year") + 
  geom_smooth(method = lm) +
  labs(title = "Max Growing Season NDVI", subtitle = "1974-2013",
       y = "NDVI\nmax(ndvi_max)") + 
  theme_bw()
max_plot


# Precipitation

# Load precipitation data for nearby site called "D1"
precip_d1 <- read.csv("data/precip_d1/d-1pdayv.ml.data.csv", na.strings = "NaN") %>% 
  mutate(DATE = as.POSIXct(date)) %>% 
  mutate(MONTH = month(DATE)) %>% 
  mutate(DAY = day(DATE)) %>% 
  mutate(YEAR = year(DATE)) %>% 
  mutate(MONTH_YEAR = paste0(MONTH, "-", YEAR)) %>% 
  mutate(HYEAR = ifelse(MONTH <= 10, YEAR, YEAR + 1)) # hydrological year starts Nov. 1 of the preceding year

str(precip_d1)

# Create new dataframe of monthly precip totals
precip_monthly_totals <- precip_d1 %>% 
  group_by(YEAR, MONTH) %>% 
  summarise(monthly_total = sum(precip, na.rm = TRUE))

# Create new dataframe of hydrological year (HY) totals
precip_HY_totals <- precip_d1 %>% 
  group_by(YEAR) %>% 
  summarize(HY_total_precip = sum(precip, na.rm = TRUE))

# Compare max NDVI with HY total precip
yearly_summaries <- merge(precip_HY_totals, landsat_ndvi_yearly_max) %>% 
  rename(year_max_ndvi = year_max)

plot_total_precip_max_ndvi <- ggplot(data = yearly_summaries, aes(x = HY_total_precip, y = year_max_ndvi)) +
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(title = "Precipitation and NDVI", 
       y = "Max NDVI for the Year", x = "Total precipitation (mm) for the hydrological year (HY)") + 
  theme_bw()
plot_total_precip_max_ndvi

# Compare spatial coefficient of variance (CV) of NDVI with HY total precip: Is there less variability in plant growth in drier or wetter years?
landsat_ndvi_yearly_min_cv <- landsat_ndvi %>% 
  group_by(YEAR) %>% 
  summarise(MIN_CV = min(ndvi_cv, na.rm = TRUE))

yearly_summaries <- merge(yearly_summaries, landsat_ndvi_yearly_min_cv)
plot_total_precip_min_cv <- ggplot(data = yearly_summaries, aes(x = HY_total_precip, y = MIN_CV)) +
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(x = "Y", y = "yearly_min_cv") + 
  theme_bw()
plot_total_precip_min_cv

# Max CV per HY as a function of precip: Is there more variability in plant growth in dryer or wetter years?
landsat_ndvi_yearly_max_cv <- landsat_ndvi %>% 
  group_by(YEAR) %>% 
  summarise(MAX_CV = max(ndvi_cv, na.rm = TRUE))

yearly_summaries <- merge(yearly_summaries, landsat_ndvi_yearly_max_cv)
plot_total_precip_max_cv <- ggplot(data = yearly_summaries, aes(x = HY_total_precip, y = MAX_CV)) +
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(title = "Spatial variability in plant growth (CV of NDVI) relative to precipitation", 
       x = "Total precipitation (mm) for the hydrological year", y = "Maximum spatial variability\nin NDVI during the growing season") + 
  theme_bw()
plot_total_precip_max_cv


# Temperature

temp_d1 <- read.csv("data/temp_d1/d-1tdayv.ml.data.csv", na.strings = "NaN") %>% 
  mutate(D1_max_temp = max_temp, na.rm = TRUE) %>% 
  mutate(D1_min_temp = min_temp, na.rm = TRUE) %>% 
  mutate(D1_mean_temp = mean_temp, na.rm = TRUE) %>% 
  mutate(date = as.POSIXct(date)) %>% 
  mutate(YEAR = year(date)) %>% 
  select("date", "D1_max_temp", "D1_min_temp", "D1_mean_temp", "YEAR")

# Summarize yearly max, min, and mean temperatures
temp_d1_yearly <- temp_d1 %>% 
  group_by(YEAR) %>% 
  summarise(max_meantemp = max(D1_mean_temp, na.rm = TRUE),
            max_maxtemp = max(D1_max_temp, na.rm = TRUE),
            min_meantemp = min(D1_mean_temp, na.rm = TRUE),
            min_mintemp = min(D1_min_temp, na.rm = TRUE),
            mean_meantemp = mean(D1_mean_temp, na.rm = TRUE))

# Calculate the growing degree days (GDD = days x temperature above freezing) for spring of each year
temp_d1_gdd_spring <- temp_d1 %>% 
  filter(month(date) < 6, D1_mean_temp > 0) %>% 
  group_by(YEAR) %>% 
  summarize(gdd_spring = sum(D1_mean_temp, na.rm = TRUE))
# Calculate the growing degree days (days x temperature above freezing) for summer of each year
temp_d1_gdd_summer <- temp_d1 %>% 
  filter(month(date) > 6, month(date) < 10, D1_mean_temp > 0) %>% 
  group_by(YEAR) %>% 
  summarize(gdd_summer = sum(D1_mean_temp, na.rm = TRUE))
temp_d1_gdd <- merge(temp_d1_gdd_spring, temp_d1_gdd_summer)  

# Merge dataframes
temp_all <- merge(temp_d1_yearly, temp_d1_gdd)
str(temp_all)

temp_all_yearly <- temp_all %>% 
  group_by(YEAR) %>% 
  summarise(TEMP_MAX_of_MEAN_Saddle_appx = (max(C1_mean_temp, na.rm = TRUE) - 10))
head(temp_all_yearly)

# Mean temp over time
ggplot(data = temp_all, aes(y = mean_meantemp, x = YEAR)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  labs(title = "Mean annual temperature over time", subtitle = "1953-2014",
       x = "Year", y = "Mean Annual Temperature (degC)") + 
  theme_bw()

# Mean spring GDD over time
ggplot(data = temp_all, aes(y = gdd_spring, x = YEAR)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  labs(title = "Spring growing degree days (GDD) over time", subtitle = "1953-2014",
       x = "Year", y = "GDD in deg C\n(baseline = 0 degC)")
temp_d1_springgdd_lm <- lm(data = temp_all, gdd_spring ~ YEAR)
summary(temp_d1_springgdd_lm)

# Mean summer GDD over time
ggplot(data = temp_all, aes(y = gdd_summer, x = YEAR)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  labs(title = "Summer growing degree days (GDD) over time", subtitle = "1953-2014",
       x = "Year", y = "GDD in deg C\n(baseline = 0 degC)")

temp_d1_summergdd_lm <- lm(data = temp_all, gdd_summer ~ YEAR)
summary(temp_d1_summergdd_lm)


# NDVI and temp
yearly_summaries2 <- merge(yearly_summaries, temp_all, by = "YEAR")
plot_meanofmean_temp_max_ndvi <- ggplot(data = yearly_summaries2, aes(x = mean_meantemp, y = year_max_ndvi)) +
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(title = "NDVI decreases with warmer mean annual temperatures",
       x = "Mean Annual Temperature (degC)", y = "Yearly Maximum NDVI")
plot_meanofmean_temp_max_ndvi