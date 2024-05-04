# project SDM to past time periods
# load packages
library(terra)
library(stringr)
library(purrr)
library(here)
library(ENMeval)
library(sf)
library(rnaturalearth)

read_stacks <- function(timeslice, bioclims = bioclims) {
  bio_ts <- rast(list.files(here("climate"), pattern = paste0("_", timeslice), full.names = TRUE))
  names(bio_ts) <- str_remove(basename(sources(bio_ts)), paste0("_", timeslice, ".tif"))
  return(bio_ts)
}



# 1. read in predictors and crop---------------------------------------------------

# read in predictors
timeslices <- paste0(0:22, "kyBP")
bioclims <- c("bio2", "bio4", "bio5", "bio8", "bio9", "bio13", "bio19")

pred_stacks <- map(timeslices, ~read_stacks(.x, bioclims = bioclims))
names(pred_stacks) <- paste0("ts_", timeslices)

# read in land boundaries for cropping
land <- rnaturalearth::ne_countries(scale = "large", continent = "South America") %>%
  vect() %>%
  terra::project("epsg:4326")

# crop and mask predictors to land boundaries
pred_stacks <- map(pred_stacks, ~crop(.x, land)) %>%
  map(~mask(.x, land))

# 2. project SDM to all time periods---------------------------------------------------

# load SDM model
sdm <- readRDS(here("output", "sdm.rds"))

# project
projections <- rast(map(pred_stacks, \(x) dismo::predict(sdm@models[["fc.LQ_rm.2"]], x)))


# 3. write out projections---------------------------------------------------
terra::writeRaster(projections, here("output", "projections.tif"), overwrite = TRUE)




