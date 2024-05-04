# script to generate the simulation time figure
library(tidyverse)
library(cowplot)
library(here)

theme_set(theme_bw())
theme_update(plot.title = element_text(size = 20),
             axis.title = element_text(size = 16),
             axis.text = element_text(size = 14),
             strip.text = element_text(size = 16),
             panel.spacing.x = unit(2, "lines"))

# 0. load the data --------------------------------------------------------
meta_files <- list.files(here("../iddc/output"), pattern = "metadata.csv", full.names = TRUE)

meta <- map(meta_files, read_csv) %>%
  map(~select(.x, -anc_sizes)) %>%
  bind_rows(.id = "file") %>%
  mutate(
    model_name = case_when(
      file == "1" ~ "One ancestral population constant size",
      file == "2" ~ "One ancestral population fluctuating size",
      file == "3" ~ "Two ancestral populations constant size",
      file == "4" ~ "Two ancestral populations fluctuating size",
    )
  )

# convert sim_time column to posix time
meta <- meta %>%
  mutate(sim_time_formatted = as.POSIXct(origin = "2024-04-04", sim_time))


mig_time <- ggplot(meta, aes(x = mig_rate, y = sim_time_formatted)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE, color = "red") +
  labs(x = "Migration rate", y = "Simulation time (Min:Sec)") +
  scale_y_datetime(date_labels = "%M:%S") +
  facet_wrap(~model_name)

n_time <- ggplot(meta, aes(x = max_local_size, y = sim_time_formatted)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE, color = "red") +
  labs(x = "Maximum local deme size", y = "Simulation time (Min:Sec)") +
  scale_y_datetime(date_labels = "%M:%S") +
  facet_wrap(~model_name)

plot_grid(mig_time, n_time, ncol = 1, labels = c("A)", "B)"))

ggsave(here("output", "plots", "sim_time_fig.png"), width = 25, height = 30, units = "cm")


