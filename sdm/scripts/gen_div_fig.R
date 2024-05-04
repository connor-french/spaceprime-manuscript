library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(patchwork)
library(here)
library(rlang)


read_gendiv <- function(folder, model) {
  paths <- list.files(here(folder), pattern = model, full.names = TRUE)

  r <- rast(paths)

  # set NA values
  r[r<=0] <- NA

  # take the mean of the stack
  r <- terra::mean(r, na.rm = TRUE)

  names(r) <- "nuc_div"

  return(r)
}


plot_gendiv <- function(gendiv, locs, title) {
  ggplot() +
    geom_spatraster(data = gendiv) +
    scale_fill_whitebox_c("muted", direction = 1, na.value = "transparent") +
    geom_sf(data = locs, shape = 21, fill = "darkgreen") +
    coord_sf(xlim = c(-55, -45), ylim = c(-31, -22)) +
    labs(title = title, fill = "Nucleotide\ndiversity") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 20),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.position = "inside",
      legend.position.inside = c(0.8, 0.3),
      legend.background = element_rect(fill = "transparent"),
      legend.key.height = unit(1.75, "lines")
    )
}

plot_pca <- function(ss_pca_df, pc_x, pc_y) {
  pc_plot <- ss_pca_df %>%
    ggplot(aes({{pc_x}}, {{pc_y}})) +
    geom_point(aes(fill = scenario), shape = 21, color = "black", size = 3) +
    scale_fill_viridis_d() +
    labs(fill = "Scenario") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 20),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 12, hjust = 0.5),
      legend.text = element_text(size = 10),
      legend.key.height = unit(2, "lines"),
      legend.position = "inside",
      legend.position.inside = c(0.8, 0.175),
      legend.background = element_rect(fill = "transparent")
    )

  return(pc_plot)
}

# 1. read in data ---------------------------------------------------------
# read in geospatial data
locs <- read_sf(here("distincta_localities_ex.geojson"))

gendiv_1p <- read_gendiv(here("../iddc/output/onepop-rasts"), "onepop_diversity_map")

gendiv_1pc <- read_gendiv(here("../iddc/output/onepop-rasts"), "onepop-map-const_diversity_map")

gendiv_2p <- read_gendiv(here("../iddc/output/twopop-rasts"), "twopop-map_diversity_map")

gendiv_2pc <- read_gendiv(here("../iddc/output/twopop-rasts"), "twopop-const-map_diversity_map")

# read in summary statistics
ss_files <- list.files(here("../iddc/output"), pattern = "sumstats.csv", full.names = TRUE)
ss <- map_df(ss_files, read_csv, .id = "file")

# PCA of sumstats
ss_pca <- prcomp(ss %>% select(-file), scale = TRUE, center = TRUE)

ss_pca_df <- as_tibble(ss_pca$x) %>%
  mutate(
    file = ss$file,
    scenario = case_when(
      file == "1" ~ "One ancestral population\nconstant size",
      file == "2" ~ "One ancestral population\nfluctuating size",
      file == "3" ~ "Two ancestral populations\nconstant size",
      file == "4" ~ "Two ancestral populations\nfluctuating size",
      .default = NA_character_
    )
  )





# 2. plot nucleotide diversity ---------------------------------------------------------
map_1p <- plot_gendiv(gendiv_1p, locs, "One ancestral population fluctuating size")

map_1pc <- plot_gendiv(gendiv_1pc, locs, "One ancestral population constant size")

map_2p <- plot_gendiv(gendiv_2p, locs, "Two ancestral populations fluctuating size")

map_2pc <- plot_gendiv(gendiv_2pc, locs, "Two ancestral populations constant size")



# 3. plot sumstat PCAs ----------------------------------------------------

pc1_pc2 <- plot_pca(ss_pca_df, pc_x = PC1, pc_y = PC2)

pc2_pc3 <- plot_pca(ss_pca_df, pc_x = PC2, pc_y = PC3) + theme(legend.position = "none")


# 4. arrange plots --------------------------------------------------------
parrange <- c(
  "
  AABB
  CCDD
  EEFF
  "
)

total_plot <- map_1pc + map_1p + map_2pc + map_2p + pc1_pc2 + pc2_pc3 +
  plot_layout(design = parrange, heights = c(1, 1)) +
  plot_annotation(tag_levels = "A", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 20))

ggsave(here("output", "plots", "gendiv_fig.png"), total_plot, width = 35, height = 50, units = "cm", dpi=300)

