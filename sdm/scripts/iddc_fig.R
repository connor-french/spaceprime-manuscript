library(tidyverse)
library(cowplot)
library(here)

theme_set(theme_bw())
theme_update(plot.title = element_text(size = 20),
             axis.title = element_text(size = 16),
             axis.text = element_text(size = 14))

# Sumstat PCA -------------------------------------------------------------

# read files
ss_files <- list.files(here("../iddc/output"), pattern = "sumstats.csv", full.names = TRUE)

ss <- map_df(ss_files, read_csv, .id = "file") %>%
  mutate(
    model_name = case_when(
      file == "1" ~ "One anc. pop. const.",
      file == "2" ~ "One anc. pop.",
      file == "3" ~ "Two anc. pop. const.",
      file == "4" ~ "Two anc. pop.",
    )
  )


# PCA of sumstats
ss_pca <- prcomp(ss %>% select(-file, -model_name), scale = TRUE, center = TRUE)

ss_pca_df <- as_tibble(ss_pca$x) %>%
  mutate(file = ss$file, model_name = ss$model_name)


# PCA plots
pc1_pc2 <- ss_pca_df %>%
  ggplot(aes(PC1, PC2, color = model_name)) +
  scale_color_viridis_d() +
  labs(color = "Model") +
  geom_point()

pc2_pc3 <- ss_pca_df %>%
  ggplot(aes(PC2, PC3, color = model_name)) +
  scale_color_viridis_d() +
  labs(color = "Model") +
  geom_point()

pc2_pc3





