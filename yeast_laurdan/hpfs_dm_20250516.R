library(here)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(ggpubr)

## Helper functions
# Function to parse and timestamp WinFLR kinetics runs
# this takes the "CSV with log" export format
parse_winflr = function(fname){
  # split on first blank line
  # file is made on Windoze w/CRLF line breaks
  sections = read_file(fname) %>% 
    str_split("\r\n\r\n", n = 2) %>% 
    unlist() %>% 
    setNames(c("fluordata", "metadata"))
  
  # parse fluorescence data
  fluordata = read_csv(
    sections[["fluordata"]], 
    skip = 2, 
    col_names = c("watch_440", "intensity_440", "watch_490", "intensity_490"),
    show_col_types = FALSE
  ) %>% 
    # convert times from min to sec
    mutate(across(contains("watch"), ~.x*60)) %>% 
    # drop columns with all NAs (caused by trailing delims)
    dplyr::select(where(~!all(is.na(.x))))
  
  # parse log, get start time
  metadata = sections[["metadata"]] %>% 
    str_split("\r\n") %>% 
    unlist() %>% 
    tibble(key = .) %>% 
    separate(key, into = c("key", "val"), sep = ':', extra = "merge", fill = "right") %>% 
    mutate(
      across(everything(), trimws),
      val = ifelse(val == "", NA, val)
    )
  
  # coregister the fluorescence readings with wall time
  fluordata %>% 
    mutate(
      watch = (watch_440 + watch_490)/2,
      clock = watch + (metadata %>% 
                         filter(key == "Collection Time") %>% 
                         distinct() %>% .$val %>% 
                         parse_date_time(orders = "%m-%d-%Y %I:%M:%S %p"))
    )
}

## Temperature calibration
calfile = here("01-hpfs", "thermalog_20250513T2028.tsv")

data_tcal = calfile %>% 
  read_tsv()

# time traces
data_tcal %>% 
  ggplot(
    aes(
      x = clock
    )
  ) +
  geom_point(aes(y = temp_int), color = "blue") +
  geom_point(aes(y = temp_act), color = "red")

# we'll call the last 5 min of every 60-min step good
data_eqbm = data_tcal %>% 
  group_by(state) %>% 
  filter(watch > min(watch) + (55*60))

# plot the cal curve data
data_eqbm %>%
  ggplot(
    aes(
      x = temp_int,
      y = temp_act
    )
  ) +
  geom_smooth(method = "lm") +
  geom_point()

# generate cal function
# note that within the range, the actual temperatures generally come out
# lower than the RTD readings bc of the lack of line loss compensation
rtd2act = function(temp_rtd){
  data_eqbm %>% 
    lm(temp_act ~ temp_int, data = .) %>%# summary()
    predict(., newdata = tibble(temp_int = temp_rtd)) %>% 
    unname()
}

# Load and collate all fluor data
files_fluor = list.files(
  here("01-hpfs"), 
  pattern = "kinetic.*.csv", 
  recursive = TRUE,
  full.names = TRUE
)

data_gp = tibble(fname = files_fluor) %>% 
  rowwise() %>%
  mutate(data = fname %>% parse_winflr() %>% list()) %>% 
  unnest(data) %>% 
  mutate(
    samp = fname %>% dirname() %>% basename(),
    gp = (intensity_440 - intensity_490)/(intensity_440 + intensity_490)
  ) %>% 
  separate(samp, into = c("lipid", "date_start"), sep = '_')

# Load and collate all T-P logs
files_tplog = list.files(
  here("01-hpfs"), 
  pattern = "tplog.*.tsv",
  full.names = TRUE
)

data_tp = files_tplog %>% read_tsv()

# generate interp functions for T and P
# we get a regularization warnings bc the clock resolution is 1 s
clock2temp = approxfun(data_tp$clock, rtd2act(data_tp$temp_int))
clock2press = approxfun(data_tp$clock, data_tp$press_act)
clock2prset = approxfun(data_tp$clock, data_tp$press_set, method = "constant", rule = 1)
clock2state = approxfun(data_tp$clock, data_tp$state, method = "constant", rule = 1)

# coregister all the data
tol_press = 1 #(1 bar)
data_gptp = data_gp %>% 
  mutate(
    temp  = clock2temp(clock),
    press = clock2press(clock),
    prset = clock2prset(clock),
    state = clock2state(clock)
  ) %>% 
  separate(lipid, into = c("cclass", "bbone"), remove = FALSE) %>% 
  mutate(
    bbone = bbone %>% factor(., levels = c(
      "APG",
      "DAG",
      "AEG",
      "DEG", 
      "W303",
      "PEM2",
      "PSD1", 
      "CHO"
    ))
  ) %>% 
  group_by(lipid, cclass, bbone) %>% 
  filter(abs(press - prset) <= tol_press) %>% 
  filter(state != max(state)) %>%
  filter(cclass %in% c("W303", "PEM2", "PSD1", "CHO")) %>%
  # keep ≤2 measurements from each state
  group_by(state, .add = TRUE) %>% 
  arrange(clock) %>% 
  filter(row_number() > n()-2) %>%
  group_by(lipid, cclass, bbone)


# run out 3D LOESS splines
res_temp  = 0.1
res_press = 1
topo_gptp = data_gptp %>%
  mutate(intensity_440 = intensity_440/max(intensity_440)) %>% 
  summarise(
    across(c(temp, press),  list(min = min, max = max)),
    mod_gp  = loess(gp ~ press + temp, span = 0.85) %>% list(),
    mod_440 = loess(intensity_440 ~ press + temp, span = 0.85) %>% list()
  ) %>% 
  rowwise() %>% 
  mutate(
    grid = expand_grid(
      temp = seq(temp_min, temp_max+res_temp, res_temp),
      press = seq(press_min, press_max+res_press, res_press)
    ) %>% list()
  ) %>% 
  unnest(grid) %>% 
  group_by(lipid, cclass, bbone) %>% 
  mutate(
    gp = predict(mod_gp[[1]], newdata = cur_data()),
    intensity_440 = predict(mod_440[[1]], newdata = cur_data())
    #se = "TRUE" #takes FOREVER for the LOESS
  )


# plot!
library(metR) # for geom_text_contour

topo_gptp %>%
  filter(press < 1000) %>%
  #filter(cclass %in% c("1bar", "250bar")) %>%
  ggplot(aes(x = temp, y = press, z = gp)) +
  facet_grid(cols = vars(cclass), scales = "free_x") +
  geom_contour_filled(aes(fill = stat(level_mid)), alpha = 0.6) +
  geom_contour(aes(color = ..level.. >= 0.5)) +
  geom_text_contour(
    stroke = 0.2,             # thin white halo for readability
    size = 2,               # text size
    label.placer = label_placer_fraction(), # keeps labels on flatter segments
    skip = 2                  # label every other contour line
  ) +
  geom_point(data = data_gptp, shape = 1, size = 0.3, alpha = 0.3) +
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +
  scale_color_manual(values = c("black", "red")) +
  theme_pubr(base_size = 6) +
  theme(
    legend.position = "right",
    legend.direction = "vertical"
  ) +
  labs(
    x = "Temperature (°C)",
    y = "Pressure (bar)",
    fill = "GP"
  )
ggsave(here("PEPCyeast_contours_gp_20250915.pdf"), width = 18, height = 6, units = "cm" )
#


topo_gpp = topo_gptp %>%
  filter(temp > 29) %>%
  filter(temp < 29.2) %>%
  filter(press < 1000) %>%
  filter(press > 0) 

data_gptp %>%
  mutate(cclass = factor(cclass, levels=c("W303", "PSD1", "CHO", "PEM2"))) %>%
  filter(temp > 28) %>%
  filter(temp < 32) %>%
  filter(press < 1000) %>%
  filter(press > 0) %>%
  ggplot(
    aes(
      x = press, 
      y = gp,
      linetype = cclass, 
      color = cclass
    )) +   
      #facet_grid(cols = vars(cclass), scales = "free_x") +
    geom_point(size = 1) +
  
    #geom_line(data = topo_gpp, size = 1) +
    scale_color_grey() +
    theme_pubclean() +
    labs(
      x = "Pressure (bar)", 
      y = "GP"
    ) +
    theme(legend.position = "top") #+
  scale_y_continuous(
    limits = c(0.0, 0.25),
    breaks = c(0.0, 0.1, 0.2)
  )
    #geom_point(size = 0.25) 

ggsave(here( "03-scratchfigs", "pressureYeast_larudanGP_dP.pdf"), width = 4, height = 4)




# bin T so replicates collapse cleanly
means_dt <- data_gptp %>%
  filter(press > 0, press < 1.5) %>%
  mutate(temp_bin = round(temp, 1)) %>%    # 0.1 °C bins; change if you want
  group_by(cclass, temp = temp_bin) %>%
  summarise(mean_gp = mean(gp, na.rm = TRUE),
            sem_gp  = sem(gp), .groups = "drop")


# slice & summarize raw data to mean ± SEM at each pressure
means_dp <- data_gptp %>%
  filter(temp > 29, temp < 29.2, press > 0, press < 1000) %>%
  mutate(press_bin = round(prset)) %>%     # or round(press, -1) for 10-bar bins
  group_by(cclass, press = press_bin) %>%
  summarise(mean_gp = mean(gp, na.rm = TRUE),
            sem_gp  = sem(gp), .groups = "drop")





topo_gpt = topo_gptp %>%
  filter(temp > 0) %>%
  #filter(temp < 29.2) %>%
  filter(press < 1.5) %>%
  filter(press > 0)

data_gptp %>%  
  mutate(cclass = factor(cclass, levels=c("W303", "PSD1", "CHO", "PEM2"))) %>%
  #c("W303", "PSD1", "CHO", "PEM2")
  #filter(temp > 30) %>%
  #filter(temp < 30.5) %>%
  filter(press > 0) %>%
  filter(press < 1) %>%
  ggplot(
    aes(
      x = temp, 
      y = gp,
      #linetype = cclass, 
      color = cclass
    )) +   
  #facet_grid(cols = vars(cclass), scales = "free_x") +
  geom_point(size = 1) +
  geom_line(data = topo_gpt, size = 1) +
  scale_color_grey() +
  theme_pubclean(base_size = 6) +
  labs(
    x = "Temperature (C)", 
    y = "C-Laurdan GP"
  ) +
  theme(legend.position = "none") #+
  scale_y_continuous(breaks=(seq(-0.3, 0.3, 0.15))) +
geom_point(size = 0.25) 

ggsave(here( "03-scratchfigs", "PressureYeast_larudanGP_dP.pdf"), width = 4, height = 4)


sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

# mean ± SEM per class at each set pressure in the 29–29.2 °C slice
gpp_sum <- data_gptp %>%
  filter(temp > 29, temp < 29.2, press > 0, press < 1000) %>%
  mutate(press_bin = round(prset)) %>%          # or: press_bin = round(press, -1)
  group_by(cclass, press = press_bin) %>%
  summarise(mean_gp = mean(gp, na.rm = TRUE),
            sem_gp  = sem(gp), .groups = "drop")

# plot: raw points (optional), LOESS fit line, then mean ± SEM
data_gptp %>%
  filter(temp > 29, temp < 29.2, press > 0, press < 1000) %>%
  ggplot(aes(x = press, y = gp, color = cclass)) +
  #geom_point(size = 0.7, alpha = 0.25) +                       # raw points, faint
  geom_line(data = topo_gpp,                                   # your fit line
            aes(x = press, y = gp, color = cclass, group = cclass),
            linewidth = 0.9) +
  geom_errorbar(data = gpp_sum,                                # mean ± SEM
                aes(x = press, ymin = mean_gp - sem_gp, ymax = mean_gp + sem_gp,
                    color = cclass),
                width = 8, linewidth = 0.35, inherit.aes = FALSE) +
  geom_point(data = gpp_sum, aes(x = press, y = mean_gp, color = cclass),
             size = 1.3, inherit.aes = FALSE) +
  scale_color_grey() +
  scale_y_continuous(limits = c(-0.1, 0.4),
  breaks = c(-0.1, 0, 0.1, 0.2, 0.3, 0.4)) +
  theme_pubclean() +
  theme(legend.position = "top") +
  labs(x = "Pressure (bar)", y = "GP")
ggsave(here( "03-scratchfigs", "PressureYeast_larudanGP_dP.pdf"), width = 2, height = 2)


data_gptp %>%
  filter(temp > 29, temp < 29.2, press > 0, press < 1000) %>%
  ggplot(aes(x = press, y = gp, color = cclass)) +
  #geom_point(size = 0.7, alpha = 0.25) +   # raw points, faint
  geom_line(data = topo_gpp,               # your fit line
            aes(x = press, y = gp, color = cclass, group = cclass),
            linewidth = 0.9) +
  geom_errorbar(data = gpp_sum,            # mean ± SEM
                aes(x = press, ymin = mean_gp - sem_gp, ymax = mean_gp + sem_gp,
                    color = cclass,),
                width = 8, linewidth = 1, inherit.aes = FALSE) +
  #geom_point(data = gpp_sum,
             #aes(x = press, y = mean_gp, color = cclass),
             #size = 1.3, inherit.aes = FALSE) +
  scale_color_grey() +
  #scale_y_continuous(limits = c(-0.3, 0.15),
                     #breaks = c(-0.3, -0.15, 0, 0.15)) +
  scale_shape_manual(values = c(17, 16, 20, 18)) +
  theme_pubclean(base_size = 6) +
  theme(legend.position = "none") +
  labs(x = "Pressure (bar)", y = "GP")


# custom greyscale values from dark to light
custom_greys <- c(
  "W303" = "#999999",   # black
  "PSD1" = "#555555",   # dark grey
  "PEM2" = "#000000",   # medium grey
  "CHO"  = "#CCCCCC"    # light grey
)

data_gptp %>%
  filter(temp > 29, temp < 30, press > 0, press < 1000) %>%
  ggplot(aes(x = press, y = gp, color = cclass)) +
  #geom_point(size = 0.7, alpha = 0.25, ) +   # raw points, faint
  geom_line(data = topo_gpp,               # fit line
            aes(x = press, y = gp, color = cclass, group = cclass),
           linewidth = 0.9) +
  geom_errorbar(data = gpp_sum,            # mean ± SEM
                aes(x = press, ymin = mean_gp - sem_gp, ymax = mean_gp + sem_gp,
                    color = cclass),
                width = 8, linewidth = 1, inherit.aes = FALSE) +
  geom_point(data = gpp_sum,
             aes(x = press, y = mean_gp, color = cclass, shape = cclass),
             size = 2, inherit.aes = FALSE) +
  scale_color_manual(values = custom_greys) +
  scale_shape_manual(values = c(0, 1, 20, 18)) +
  scale_y_continuous(limits = c(-0.3, 0.15),
                     breaks = c(-0.3, -0.15, 0, 0.15)) +
  theme_pubclean() +
  theme(legend.position = "top") +
  labs(x = "Pressure (bar)", y = "GP")
ggsave(here( "03-scratchfigs", "PEPCYeast_larudanGP_dP.pdf"), width = 4, height = 4)

sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

# custom greys (match pressure plothttp://127.0.0.1:40373/graphics/plot_zoom_png?width=333&height=259)
custom_greys <- c(
  "W303" = "#000000",   # black
  "PSD1" = "#555555",   # dark grey
  "PEM2" = "#999999",   # medium grey
  "CHO"  = "#CCCCCC"    # light grey
)

# mean ± SEM for temperature sweep at ~1 bar
means_dt <- data_gptp %>%
  filter(press > 0, press <1) %>%
  mutate(
    press = 1, 
    temp_bin = round(temp, 20)) %>%  # 1 °C bins
  group_by(cclass, temp = temp_bin) %>%
  summarise(mean_gp = mean(gp, na.rm = TRUE),
            sem_gp  = sem(gp), .groups = "drop")

# LOESS fit subset for ~1 bar
fits_dt <- topo_gpt   # from your earlier code

# plot
data_gptp %>%
  filter(press > 0, press < 1) %>%
  mutate(press = 1) %>%
  ggplot(aes(x = temp, y = gp, color = cclass)) +
  #geom_point(size = 0.7, alpha = 0.25) +                                   # raw points
  #geom_line(data = fits_dt,
   #         aes(x = temp, y = gp, color = cclass, group = cclass),
    #        linewidth = 0.9) +                                             # fit line
  geom_pointrange(data = means_dt,
                aes(x = temp, ymin = mean_gp - sem_gp, ymax = mean_gp + sem_gp,
                    color = cclass, size = 5)) +        # thick black error bars
  #geom_point(data = means_dt,
            # aes(x = temp, y = mean_gp, color = cclass, shape = cclass),
             #size = 2, inherit.aes = FALSE) +
  scale_shape_manual(values = c(17, 16, 20, 18)) +
  scale_color_manual(values = custom_greys) +
  scale_y_continuous(limits = c(-0.3, 0.15),
                     breaks = c(-0.3, -0.15, 0, 0.15)) +                   # same scale as pressure plot
  theme_pubclean() +
  theme(legend.position = "top") +
  labs(x = "Temperature (°C)", y = "C-Laurdan GP")
ggsave(here( "03-scratchfigs", "PEPCYeast_larudanGP_dT.pdf"), width = 4, height = 4)

# choose reference
ref_class <- "W303"

# join fits with SEM
diff_df <- fits_dt_se %>%            # your LOESS fit ± SEM data frame
  select(cclass, temp, gp_fit, gp_sem) %>%
  inner_join(
    fits_dt_se %>% filter(cclass == ref_class) %>%
      select(temp, gp_fit_ref = gp_fit, gp_sem_ref = gp_sem),
    by = "temp"
  ) %>%
  filter(cclass != ref_class) %>%
  mutate(
    diff   = gp_fit - gp_fit_ref,
    seDiff = sqrt(gp_sem^2 + gp_sem_ref^2),
    lo95   = diff - 1.96 * seDiff,
    hi95   = diff + 1.96 * seDiff
  )

# plot the difference
ggplot(diff_df, aes(x = temp, y = diff, fill = cclass, color = cclass)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.2, color = NA) +
  geom_line() +
  labs(y = paste0("ΔGP vs ", ref_class, " (fit ± 95% CI)"),
       x = "Temperature (°C)") +
  theme_pubclean()


# plot: raw points (optional), LOESS fit line, then mean ± SEM
data_gptp %>%
  filter(temp > 29, temp < 29.2, press > 0, press < 1000) %>%
  ggplot(aes(x = press, y = gp, color = cclass, shape = cclass)) +
  geom_line(
    data = topo_gpp,
    aes(x = press, y = gp, color = cclass, group = cclass),
    linewidth = 1
  ) +
  geom_errorbar(data = gpp_sum,            # mean ± SEM
                aes(x = press, ymin = mean_gp - sem_gp, ymax = mean_gp + sem_gp,
                    color = cclass), width = 4, linewidth = 1, inherit.aes = FALSE) +
  geom_point(data = gpp_sum,
             aes(x = press, y = mean_gp, color = cclass, shape = cclass),
             size = 3, inherit.aes = FALSE) +
  scale_color_manual(values = custom_greys) +
  scale_shape_manual(values = c(17, 16, 20, 18)) +
  scale_y_continuous(limits = c(-0.3, 0.15),
                     breaks = c(-0.3, -0.15, 0, 0.15)) +# adjust count/order to your cclass levels
  theme_pubclean() +
  theme(legend.position = "none") +
  labs(x = "Pressure (bar)", y = "GP")
ggsave(here( "03-scratchfigs", "newPressureYeast_larudanGP_dP.pdf"), width = 3, height = 3)

  data_gptp %>%
  filter(temp > 29, temp < 29.2, press > 0, press < 1000) %>%
  ggplot(aes(x = press, y = gp, color = cclass, shape = cclass)) +
  geom_point(size = 0.7, alpha = 0.25) +   # raw points, faint
  geom_line(data = topo_gpp,               # your fit line
            aes(x = press, y = gp, color = cclass, group = cclass),
            linewidth = 0.9, show.legend = TRUE) +
  geom_errorbar(data = gpp_sum,            # mean ± SEM
                inherit.aes = FALSE,
                aes(x = press, ymin = mean_gp - sem_gp, ymax = mean_gp + sem_gp),
                color = cclass, width = 8, linewidth = 1) +    # thick black bars
  geom_point(data = gpp_sum,
             aes(x = press, y = mean_gp, color = cclass, shape = cclass),
             size = 2, inherit.aes = FALSE) +
  scale_color_grey() +
  scale_shape_manual(values = c(17, 16, 20, 18)) +
  scale_y_continuous(limits = c(-0.3, 0.15),
                     breaks = c(-0.3, -0.15, 0, 0.15)) +
  theme_pubclean() +
  theme(legend.position = "top") +
  labs(x = "Pressure (bar)", y = "GP")


# choose reference
ref_class <- "W303"

# join fits with SEM
diff_df <- fits_dt_se %>%            # your LOESS fit ± SEM data frame
  select(cclass, temp, gp_fit, gp_sem) %>%
  inner_join(
    fits_dt_se %>% filter(cclass == ref_class) %>%
      select(temp, gp_fit_ref = gp_fit, gp_sem_ref = gp_sem),
    by = "temp"
  ) %>%
  filter(cclass != ref_class) %>%
  mutate(
    diff   = gp_fit - gp_fit_ref,
    seDiff = sqrt(gp_sem^2 + gp_sem_ref^2),
    lo95   = diff - 1.96 * seDiff,
    hi95   = diff + 1.96 * seDiff
  )
sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

# helper for LOESS fit with SEM enabled
fit_loess <- function(df, span = 0.4) {
  loess(gp ~ press + temp,
        data = df,
        span = span,
        control = loess.control(surface = "interpolate"))
}

predict_with_se <- function(mod, newdata) {
  pr <- predict(mod, newdata = newdata, se = TRUE)
  tibble(
    gp_fit = as.numeric(pr$fit),
    gp_sem = as.numeric(if (!is.null(pr$se.fit)) pr$se.fit else pr$se)
  )
}

# create a grid for prediction at ~1 bar
grid_dt <- expand.grid(
  temp  = seq(min(data_gptp$temp, na.rm = TRUE),
              max(data_gptp$temp, na.rm = TRUE),
              length.out = 200),
  press = 1
)

# fit LOESS for each cclass and predict ± SEM
fits_dt_se <- data_gptp %>%
  filter(press > 0, press < 1.5) %>%
  group_by(cclass) %>%
  group_modify(~{
    mod <- fit_loess(.x)
    bind_cols(cclass = .y$cclass, grid_dt, predict_with_se(mod, grid_dt))
  }) %>%
  ungroup()

# plot the difference
ggplot(diff_df, aes(x = temp, y = diff, fill = cclass, color = cclass)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.2, color = NA) +
  geom_line() +
  labs(y = paste0("ΔGP vs ", ref_class, " (fit ± 95% CI)"),
       x = "Temperature (°C)") +
  theme_pubclean()



ggsave(here( "03-scratchfigs", "PEPCYeast_larudanGP_dP.pdf"), width = 4, height = 4)

sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

# ---- summarize means ± SEM vs Temperature at ~1 bar ----
means_dt <- data_gptp %>%
  filter(press > 0, press < 1.5) %>%
  mutate(temp_bin = round(temp / 5) * 5) %>%   # rounds to nearest 5 °C
  group_by(cclass, temp = temp_bin) %>%
  summarise(mean_gp = mean(gp, na.rm = TRUE),
            sem_gp  = sem(gp), .groups = "drop")

# fitted slice for ~1 bar (use your existing fit grid)
fits_dt <- topo_gpt   # from your earlier code (press in (0, 1.5))
# If you built topo_with_se with gp_fit/gp_sem, you can swap it in below.

# ---- plot: raw points (faint), LOESS fit line, then mean ± SEM ----
data_gptp %>%
  filter(press > 0, press < 1.5) %>%
  #filter(temp > 29, temp < 29.2, press > 0, press < 1000) %>%
  ggplot(aes(x = temp, y = gp, color = cclass)) +
  #geom_point(size = 0.7, alpha = 0.25) +                                   # raw points
  geom_line(data = fits_dt,
            aes(x = temp, y = gp, color = cclass, group = cclass),
            linewidth = 1, show.legend = TRUE) +                          # LOESS fit
  geom_errorbar(data = means_dt,
                aes(x = temp, ymin = mean_gp - sem_gp, ymax = mean_gp + sem_gp,
                    color = cclass),
                width = 1, linewidth = 1, inherit.aes = FALSE) +      # mean ± SEM
  geom_point(data = means_dt,
             aes(x = temp, y = mean_gp, color = cclass, shape = cclass),
             size = 3, inherit.aes = FALSE) +
  scale_color_manual(values = custom_greys) +
  scale_shape_manual(values = c(17, 16, 20, 18)) +
  theme_pubclean() +
  theme(legend.position = "none") +
  labs(x = "Temperature (°C)", y = "C-Laurdan GP") +
  scale_y_continuous(limits = c(-0.3, 0.15),
                   breaks = c(-0.3, -0.15, 0, 0.15))
ggsave(here( "03-scratchfigs", "newPressureYeast_larudanGP_dT.pdf"), width = 3, height = 3)

