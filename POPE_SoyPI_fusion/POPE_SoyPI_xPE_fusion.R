# ---- Packages ----
library(tidyverse)
library(broom)

# ---- Load data ----
df <- read_csv("C:/Users/danie/OneDrive/Desktop/Pressure/POPE_SoyPI_fusion/fusion_POPE_SoyPI.csv"
)

# ---- Tidy to long format: composition + replicate ----
df_long <- df %>%
  pivot_longer(
    cols = -`time (s)`,
    names_to   = c("composition","replicate"),
    names_pattern = "([A-Za-z0-9]+)\\s*([ABC])"
  ) %>%
  rename(time = `time (s)`, intensity = value) %>%
  mutate(
    composition = factor(composition, levels = c("POPE","SoyPI","xPE", "POPG"))
  )

# ---- Mean ± SEM at each time for each composition ----
sum_dt <- df_long %>%
  group_by(composition, time) %>%
  summarise(
    n   = sum(!is.na(intensity)),
    mean_int = mean(intensity, na.rm = TRUE),
    sd_int   = sd(intensity,   na.rm = TRUE),
    sem_int  = sd_int / sqrt(n),
    .groups = "drop"
  )

# ---- Single-exponential fits to the mean curves ----
# Model: y = A * (1 - exp(-k * t))
fits <- sum_dt %>%
  group_by(composition) %>%
  group_modify(~{
    dat <- .x %>% filter(is.finite(time), is.finite(mean_int))
    if (nrow(dat) < 5) return(tibble())
    
    A0 <- max(dat$mean_int, na.rm = TRUE)
    t_half <- dat$time[which.min(abs(dat$mean_int - A0/2))]
    k0 <- if (is.finite(t_half) && t_half > 0) log(2)/t_half else 0.01
    
    fit <- try(
      nls(mean_int ~ A * (1 - exp(-k * time)),
          data = dat,
          start = list(A = A0, k = k0),
          control = nls.control(maxiter = 500, warnOnly = TRUE)),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) return(tibble())
    
    rng <- range(dat$time, finite = TRUE)
    newx <- tibble(time = seq(rng[1], rng[2], length.out = 400))
    
    augment(fit, newdata = newx) %>%
      transmute(time, fit_int = .fitted)
  }) %>%
  ungroup() %>%
  mutate(composition = rep(levels(df_long$composition), each = 400)[1:n()])

# ---- Plot: all compositions in a single graph ----
ggplot(sum_dt, aes(x = time, y = mean_int*(100), color = composition, fill = composition)) +
  geom_ribbon(aes(ymin = (mean_int - sem_int)*100, ymax = (mean_int + sem_int)*100),
              alpha = 0.30, color = NA) +
  geom_line(linewidth = 0.5) +
  geom_line(data = fits, aes(y = fit_int*100), linewidth = 0.9, linetype = 2, show.legend = FALSE) +
  labs(
    #title = "FRET Liposomal Fusion (Mean ± SEM)",
    x = "Time (s)",
    y = "Fusion efficiency (%)",
    color = "Composition",
    fill  = "Composition"
  ) +
  theme_classic(base_size = 6) +
  theme(legend.position = "none")


ggsave("C:/Users/danie/OneDrive/Desktop/Pressure/POPE_SoyPI_fusion/POPE_SoyPI_xxPE_POPG.pdf", width = 4, height = 4, units = 'cm', dpi = 300)  # Save as PDF

