library(tidyverse)
library(terra)
library(sf)
library(tidyterra)
library(ggspatial)
library(cowplot)
library(rnaturalearth)
library(here)


# functions ---------------------------------------------------------------

sum_projection <- function(pnum, proj = projections, max_local_size = 1000) {
  proj1 <- proj[[pnum]]
  sum <- sum(values(proj1), na.rm = TRUE) * max_local_size
  return(sum)
}


# 0. load the data --------------------------------------------------------

projections <- rast(here("output", "projections.tif"))

locs <- read_sf(here("distincta_localities_ex.geojson"))

af <- read_sf(here("climate", "atlantic_forest.geojson"))

# range outline for inset plot
outline <- projections$ts_0kyBP %>%
  as.polygons() %>%
  aggregate()

# 2. calculate change through time ----------------------------------------

size_df <- map_dbl(1:nlyr(projections), sum_projection) %>%
  as_tibble_col(column_name = "size") %>%
  mutate(year = 0:22)




# 3. plot ----------------------------------------------------------------

# size change through time plot
line_plot <- ggplot(size_df, aes(x = year, y = size)) +
  geom_line() +
  geom_point() +
  labs(x = "Years before present (kya)", y = "Total landscape size") +
  scale_y_continuous(labels = scales::scientific) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

## inset map
brazil <- ne_countries(scale = "medium", country = "brazil", returnclass = "sf")

brazil_inset <- ggplot() +
  geom_sf(data = brazil, fill = "white") +
  geom_sf_text(data = brazil, aes(label = name), size = 5.5, nudge_x = -2.55, nudge_y = 4.25) +
  #geom_spatvector(data = af, fill = "transparent", color = "darkgreen") +
  geom_spatvector(data = outline, fill = "darkgreen", color = "darkgreen") +
  theme_void()


map_0 <- ggplot() +
  geom_spatraster(data = projections$ts_0kyBP) +
  scale_fill_whitebox_c("muted", limits = c(0, 1)) +
  geom_sf(data = locs, shape = 21, fill = "darkgreen") +
  coord_sf(xlim = c(-55, -45), ylim = c(-31, -22)) +
  # add a north arrow
  annotation_north_arrow(which_north = "true") +
  annotation_scale(location = "br", width_hint = 0.35) +
  labs(title = "Present") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  ) +
  annotation_custom(
    grob = ggplotGrob(brazil_inset),
    xmin = -48.5, xmax = -44.9, ymin = -31.5, ymax = -27
  )


map_6 <- ggplot() +
  geom_spatraster(data = projections$ts_6kyBP) +
  scale_fill_whitebox_c("muted", limits = c(0, 1)) +
  coord_sf(xlim = c(-55, -45), ylim = c(-31, -22)) +
  labs(title = "Mid-Holocene (6 kya)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

map_22 <- ggplot() +
  geom_spatraster(data = projections$ts_22kyBP) +
  scale_fill_whitebox_c("muted", limits = c(0, 1)) +
  coord_sf(xlim = c(-55, -45), ylim = c(-31, -22)) +
  labs(
    title = "Last Glacial Maximum (22 kya)",
    fill = "Habitat\nSuitability"
    ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.82, 0.25),
    legend.background = element_blank(),
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )



# Set up the plot grid
top_row <- plot_grid(map_0, map_6, map_22, ncol = 3)
bottom_row <- plot_grid(line_plot, ncol = 1)

total_plot <- plot_grid(
  top_row,
  bottom_row,
  nrow = 2,
  rel_heights = c(1, 0.7),
  rel_widths = c(1, 0.6)
)

# save plot
ggsave(here("output", "plots", "sdm-fig.svg"), total_plot, width = 40, height = 20, units = "cm")




