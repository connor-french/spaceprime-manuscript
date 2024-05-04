
# List of packages to be installed
packages <- c("terra", "dplyr", "stringr", "sf", "here", "readr")

# Function to check and install packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Apply the function to each package
sapply(packages, install_if_missing)


library(terra)
library(dplyr)
library(stringr)
library(sf)
library(here)
library(readr)

pull_climate <- function(in_file, out_file, agg, func="mean", crop_mcp) {
  r <- terra::rast(in_file)

  v <- terra::aggregate(r, fact=agg, fun=func)

  c <- terra::crop(v, crop_mcp) %>%
    terra::mask(crop_mcp)

  writeRaster(c, out_file, overwrite=TRUE)

}


download_traceclim <- function(locs_sf, buffer, agg, bioclims, outpath) {

  mcp_all <- st_convex_hull(st_union(locs_sf)) %>%
    st_buffer(dist = units::set_units(buffer, degree)) %>%
    terra::vect()

  time_id <- seq(-200, 10, 10)
  ky_bp <- paste0(seq(22, 1, -1), "kyBP")

  for (i in bioclims) {
    for (j in 1:length(time_id)) {
      success <- FALSE
      attempts <- 0
      while (!success && attempts < 3) {
        attempts <- attempts + 1
        bioclim_file <- paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/chelsa_trace/bio/CHELSA_TraCE21k_", i, "_", time_id[[j]], "_V1.0.tif")
        bioclim_out <- file.path(outpath, paste0(i, "_", ky_bp[[j]], ".tif"))

        tryCatch({
          pull_climate(in_file = bioclim_file, out_file = bioclim_out, agg = agg, crop_mcp = mcp_all)
          success <- TRUE
        }, error = function(e) {
          message(sprintf("Attempt %d failed: %s", attempts, e$message))
          if (attempts == 3) {
            message(sprintf("Skipping download after %d attempts for %s", attempts, bioclim_file))
          }
        })
      }
    }
  }
}

download_current <- function(locs_sf, buffer, agg, bioclims, outpath) {

  mcp_all <- st_convex_hull(st_union(locs_sf)) %>%
    st_buffer(dist = units::set_units(buffer, degree)) %>%
    terra::vect()

  for (i in bioclims) {
    success <- FALSE
    attempts <- 0
    while (!success && attempts < 3) {
      attempts <- attempts + 1

      bioclim_file <- paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_", i, "_1981-2010_V.2.1.tif")
      bioclim_out <- file.path(outpath, paste0(i, "_0kyBP.tif"))

      tryCatch({
        pull_climate(in_file = bioclim_file, out_file = bioclim_out, agg = agg, crop_mcp = mcp_all)
        success <- TRUE
      }, error = function(e) {
        message(sprintf("Attempt %d failed: %s", attempts, e$message))
        if (attempts == 3) {
          message(sprintf("Skipping download after %d attempts for %s", attempts, bioclim_file))
        }
      })
    }
  }
}

bioclims_curr <- c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9")

bioclims_trace <- c("bio02", "bio04", "bio05", "bio08", "bio09", "bio13", "bio19")

locs_sf <- read_csv(
  here("distincta_localities.csv")
) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

#download_current(locs_sf = locs_sf, buffer = 2, agg = 10, bioclims = bioclims_curr, outpath = here("climate"))
download_traceclim(locs_sf = locs_sf, buffer = 2, agg = 10, bioclims = bioclims_trace, outpath = here("climate"))


