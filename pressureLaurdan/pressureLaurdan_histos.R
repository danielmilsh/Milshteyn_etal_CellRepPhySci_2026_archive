library(tidyverse)
library(readr)

# ---- 1) Read data
df <- read_csv("~/pressureLaurdan.csv", show_col_types = FALSE)

# If your first column name is exactly "GP values" (as in your file):
df_long <- df %>%
  rename(gp = `GP values`) %>%
  pivot_longer(
    cols = -gp,
    names_to = "sample",
    values_to = "density"   # this is whatever your histogram y-values represent
  ) %>%
  separate(sample, into = c("cell", "pressure"), sep = "_") %>%
  mutate(
    pressure = factor(pressure, levels = c("1", "250"),
                      labels = c("1 bar", "250 bar")),
    cell = factor(cell, levels = c("W303", "FM628", "K562"))
  )

# ---- 2) Plot: overlay curves, facet by cell type
ggplot(df_long, aes(x = gp, y = density, color = pressure)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ cell, ncol = 1, scales = "fixed") +
  labs(
    x = "Laurdan GP",
    y = "Histogram (density / fraction)",
    color = "Pressure"
  ) +
  theme_classic(base_size = 6) +
  scale_color_grey(start = 0.2, end = 0.7) +
  scale_color_manual(values = c(
    "1 bar"   = "#08306B",   # dark blue
    "250 bar" = "#7F0000"    # dark red
  )) +
  theme(legend.position = "none") 

