# sampling table
library(geobr)
library(tidyverse)
library(gt)
library(sf)
library(here)
library(glue)



# 0. read in data ---------------------------------------------------------

# read in data
locs <- read_sf(here("distincta_localities_ex.geojson"))

# download municipalities
mun <- read_municipality(code_muni = "all", year = 2020)

# extract municipality information
locs_mun <-
  locs %>%
  # transform to same crs as municipalities
  st_transform(st_crs(mun)) %>%
  st_join(mun, join = st_within) %>%
  # transform back to 4326
  st_transform(4326) %>%
  # translate ancestral population id to text
  mutate(
    anc_pop_id = glue("ANC{anc_pop_id}")
  )


# 1. create sampling table -----------------------------------------------
# create sampling table
samp_table_supp <- locs_mun %>%
  as_tibble() %>%
  select(
    "Individual ID" = individual_id,
    "Ancestral Population" = anc_pop_id,
    "Latitude" = latitude,
    "Longitude" = longitude,
    "Municipality" = name_muni,
    "State" = name_state
  ) %>%
  # arrange by ancestral population ID, state, then municipality
  arrange(`Ancestral Population`, State, Municipality) %>%
  gt() %>%
  fmt_number(
    columns = c(Latitude, Longitude),
    decimals = 2
  ) %>%
  # add a footnote for latitude and longitude
  tab_footnote(
    footnote = "Latitude and longitude are in decimal degrees.",
    locations = cells_column_labels(columns = c(Latitude, Longitude))
  ) %>%
  cols_align(
    align = "center",
    columns = c(Latitude, Longitude, Municipality, State)
  ) %>%
  cols_width(
    c(Latitude, Longitude, Municipality, State) ~ px(100)
  )

gtsave(samp_table_supp, here("output", "tables", "supp_sampling_table.tex"))
gtsave(samp_table_supp, here("output", "tables", "supp_sampling_table.docx"))



