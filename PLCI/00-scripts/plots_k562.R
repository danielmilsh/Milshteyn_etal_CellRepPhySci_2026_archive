library(tidyverse)
library(here)
library(gridExtra)
library(httr)

source(here("00-scripts", "seelipids_helpers.R"))
source(here("00-scripts", "parse_lipidmaps.R"))
source(here("00-scripts", "c0op_helpers.R"))

# where to store metadata
file_metdat = here("01-metadata", "metadata.csv")

# load the parsed data
file_tidy = here("02-tidydata", "lipidmaps_tidy.tsv")
# skip if parsing was just done and lmapdata_long is already in the global environment
lmapdata_abs = read_tsv(file_tidy)

# this is a Google Sheets key for the metadata file, to facilitate a nice Excel-like interface
gsht_metdat = "1w06ib3W2fOUGj7iOjV9cvEZ4mWMnG1YH1deafM71KgQ"

# refresh c0 spreadsheet from online interface
  GET(str_glue("https://docs.google.com/spreadsheet/ccc?key={gsht_metdat}&output=csv")) %>% 
  content("raw") %>% 
  write_file(file_metdat)

# and load the data in
metadata = file_metdat %>% read_csv()

# join metadata to lipid data
lmapdata_raw = lmapdata_abs %>%
  left_join(metadata, by = "eid") %>% 
  # then make compound ID a factor so it plots in a uniform order
  mutate(
    # put headgroups in the order specified in the color palette,
    # which can be changed in seelipids_helpers.R
    class = factor(class, levels = names(chroma_cl)),
    # then order things by total chain length and unsaturation
    id = factor(
      id,
      levels = cur_data() %>%
        distinct(id, class, carbon, dbonds) %>%
        arrange(class, carbon, dbonds) %>%
        pull(id) %>%
        unique()
    )
  ) %>% 
  # sort the rows for easy viewing
  arrange(pres_gro, time_gro, rep) %>% 
  mutate(eid = eid %>% factor(levels = unique(eid))) %>% 
  arrange(eid, id)

## we have 39:0 and 39:1 PGs at fairly high abundance
## but I am gonna leave them alone for now
#unfragd = lmapdata_raw %>% 
#  filter(id == "PE 41:6_19.26") %>% 
#  select(gn, sp, eid, class, frac_molar, dbonds, carbon, temp_stor) %>% 
#  cross_join(
#    tibble(
#      id = "PE 19:0/22:6_19.26",
#      annot = "19:0/22:6",
#      carbsn1 = 19,
#      carbsn2 = 22,
#      dbonsn1 = 0,
#      dbonsn2 = 6
#    )
#  )

# other than that we are only considering species with MS2-confirmed stereochem.
lmapdata = lmapdata_raw %>% 
  filter(ms2)# %>% 
  #bind_rows(unfragd)

# store it
lmapdata %>% write_tsv(here("02-tidydata", "lipidmaps_wmeta.tsv"))

#### PLOTS

pldata = lmapdata %>% 
  arrange(eid) %>% 
  # here's where we filter to phospholipids only
  group_by(eid) %>%
  filter(str_detect(class, 'P') & !str_detect(class, 'Cer')) %>% 
  filter(class %in% c("LPC", "LPE", "PS", "PC", "O-PC", "P-PC", "PI", "PE", "P-PE", "O-PE")) %>%
  # and renormalize to total phospholipids
  # this relies on the above group_by() call
  mutate(frac_molar = frac_molar / sum(frac_molar))

pldata %>% ungroup() %>% distinct(class, annot) %>% nrow()


pldata %>% 
  #filter(class %in% c("P-PE", "O-PE")) %>%
  gg_headgp(
    darkmode = FALSE,
    # aesthetic mappings are passed straight thru
    x = eid,
    y = frac_molar,
    fill = class
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(
    title = "K562 phospholipids",
    x = "Sample ID",
    y = "Mole fraction phospholipids",
    fill = "Lipid class"
  )
# save vector and raster images
ggsave(here("04-pdf", "k562_phospholipids_ms2_20250618a.pdf"), width = 8, height = 5)
ggsave(here("05-png", "k562_phospholipids_ms2_20250618a.png"), width = 8, height = 5)

class_foci = c("PC", "O-PC", "P-PC" , "PE", "P-PE", "O-PE", "PI", "PS")

# acyl chain comparisons
pldata %>% 
  filter(class %in% class_foci) %>% 
  arrange(eid) %>% 
  # now add all the lipids in every class and #C together
  group_by(sp, eid, pres_gro, time_gro, class, carbon) %>% 
  summarize(frac_molar = sum(frac_molar)) %>% 
  # now average by genus
  group_by(sp, pres_gro, time_gro, class, carbon) %>% 
  summarize(frac_molar = mean(frac_molar)) %>%
  ungroup() %>% 
  arrange(pres_gro, time_gro) %>% 
  mutate(
    trt = str_glue("{pres_gro} bar, {time_gro} h"),
    trt = trt %>% factor(levels = unique(.))
  ) %>% 
  gg_acylch(
    darkmode = FALSE,
    meanline = TRUE,
    # this shades the odd columns for contrast
    alpha_odd = 0.5,
    x = carbon,
    y = frac_molar,
    fill = class,
    rows = vars(trt)
  ) +
  #coord_flip() +
  labs(
    title = "K562 chain lengths",
    x = "Total acyl carbons",
    y = "Relative frequency"
  )
# save vector and raster images
ggsave(here("04-pdf", "k562_ms2_acylcarbons_20250625a.pdf"), width = 8, height = 4)
ggsave(here("05-png", "k562_ms2_acylcarbons_20250625a.png"), width = 8, height = 4)

pldata %>% 
  filter(class %in% class_foci) %>% 
  arrange(eid) %>% 
  # now add all the lipids in every class and #C together
  group_by(sp, eid, pres_gro, time_gro, class, dbonds) %>% 
  summarize(frac_molar = sum(frac_molar)) %>% 
  # now average by genus
  group_by(sp, pres_gro, time_gro, class, dbonds) %>% 
  summarize(frac_molar = mean(frac_molar)) %>%
  # normalize within each headgroup
  ungroup() %>% 
  arrange(pres_gro, time_gro) %>% 
  mutate(
    trt = str_glue("{pres_gro} bar, {time_gro} h"),
    trt = trt %>% factor(levels = unique(.))
  ) %>% 
  gg_acylch(
    darkmode = FALSE,
    meanline = TRUE,
    # this shades the odd columns for contrast
    alpha_odd = 1,
    x = dbonds,
    y = frac_molar,
    fill = class,
    rows = vars(trt)
  ) +
  #coord_flip() +
  labs(
    title = "K562 unsaturations",
    x = "Total unsaturations",
    y = "Relative frequency"
  )
# save vector and raster images
ggsave(here("04-pdf", "k562_ms2_acyldbonds_20250625a.pdf"), width = 8, height = 4)
ggsave(here("05-png", "k562_ms2_acyldbonds_20250625a.png"), width = 8, height = 4)

# Thin out the PL data by removing any lipid class that is <0.5 mol% in all samples
class_minor = pldata %>% 
  group_by(class, eid) %>% 
  summarise(frac_molar = sum(frac_molar)) %>% 
  filter(max(frac_molar) < 0.00005) %>% 
  .$class %>% unique() %>% print()

# Estimate curvature contributions
# runs all 3 schemes by default
plcidata = pldata %>% 
  filter(!(class %in% class_minor)) %>% 
  calc_plci() %>%
  filter(scheme == "linreg")

CORR = -0.00875896

library(dplyr)
library(stringr)

# ---- SETTINGS (yours) ----
CORR <- -0.00875896

CURV_COL   <- "c0"          # per-lipid curvature
ABUND_COL  <- "frac_molar"  # mol fraction (0–1)
SN2_DB_COL <- "dbonsn2"     # number of DB on sn2 (or we parse from 'sn2')
GROUP_VARS <- c("eid", "id")

# ---- HELPER ----
# If 'dbonsn2' is missing, try to derive from a string column 'sn2' like "22:6".
derive_sn2_dbonds <- function(df, sn2_db_col = SN2_DB_COL) {
  if (!sn2_db_col %in% names(df)) {
    if ("sn2" %in% names(df)) {
      df <- df %>%
        mutate(sn2_dbonds = as.numeric(str_extract(sn2, "(?<=:)\\d+")))
    } else {
      stop("Neither 'dbonsn2' nor a parsable 'sn2' (e.g., '22:6') is present.")
    }
  } else {
    df <- df %>%
      mutate(sn2_dbonds = as.numeric(.data[[sn2_db_col]]))
  }
  df
}

# ---- WORKFLOW ----
plcidata2 <- plcidata %>%
  # ensure we have a numeric sn2_dbonds column
  derive_sn2_dbonds() %>%
  # flag PUFA with >=2 DB on sn2
  mutate(sn2_is_pufa2p = !is.na(sn2_dbonds) & sn2_dbonds >= 2) %>%
  # corrected per-lipid curvature
  mutate(
    curvature_corr = as.numeric(.data[[CURV_COL]]) + ifelse(sn2_is_pufa2p, CORR, 0),
    frac_molar     = as.numeric(.data[[ABUND_COL]]),
    # SAFE product (avoid NA poisoning)
    plci      = coalesce(curvature_corr, 0) * coalesce(frac_molar, 0)
  )


# plot PLCI by individual
plcidata2 %>% 
  filter(scheme == "linreg") %>% 
  arrange(eid) %>% 
  group_by(eid, class, scheme) %>% #filter(class=="LPC") %>% View()
  #filter(plci<0) %>% 
  summarise(plci = sum(plci)) %>% 
  gg_plcurv(
    darkmode = FALSE,
    #rows = scheme,
    cols = eid,
    y = plci
  ) +
  # make negative go upward!
  scale_y_reverse() +
  theme_pubr() +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    strip.text.x = element_text(angle = 90),
    legend.position = "right"
  ) +
  labs(
    title = "K562 PL curvature by sample",
    fill = "Lipid class",
    x = "Sample"
  )

ggsave(here("04-pdf", "k562_ms2_plcisamp_20250626a.pdf"), width = 6, height = 4)
ggsave(here("05-png", "k562_ms2_plcisamp_20250626a.png"), width = 6, height = 4)

# plot PLCI mean/SEM
plcidata2 %>% 
  filter(scheme == "linreg") %>% 
  filter(time_gro == 6) %>%
  arrange(eid) %>% 
  group_by(eid, pres_gro, time_gro, class, scheme) %>% 
  summarise(plci = sum(plci)) %>% 
  # average by treatment and also get SEMs
  group_by(eid, pres_gro, time_gro, scheme) %>%
  mutate(plci_tot = ifelse(row_number() == 1, sum(plci), NA)) %>% 
  group_by(pres_gro, time_gro, class, scheme) %>%
  summarize(
    plci = mean(plci),
    plci_sem = sd(plci_tot, na.rm = TRUE) / sqrt(n()),
    plci_tot = mean(plci_tot, na.rm = TRUE)
  ) %>% 
  arrange(pres_gro, time_gro) %>% 
  mutate(
    trt = str_glue("{pres_gro} bar, {time_gro} h"),
    trt = trt %>% factor(levels = unique(.))
  ) %>% 
  gg_plcurv(
    darkmode = FALSE,
    #rows = scheme,
    cols = trt,
    y = plci
  ) +
  # there are two errorbar colors
  # so you can make the downward-pointing side white-on-black if you like
  geom_errorbar(
    aes(
      x = 0.3, # hardcoded into gg_plcurv()
      ymin = plci_tot - plci_sem,
      ymax = plci_tot
    ),
    color = "black",
    #color = "white",
    width = 0
  ) +
  geom_errorbar(
    aes(
      x = 0.3, # hardcoded into gg_plcurv()
      ymin = plci_tot,
      ymax = plci_tot + plci_sem
    ),
    color = "black",
    width = 0
  ) +
  # make negative go upward!
  #scale_y_reverse() +
  theme_pubr() +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    strip.text.x = element_text(angle = 90),
    legend.position = "none"
  ) +
  labs(
    fill = "Lipid class",
    x = "Treatment"
  ) +
  scale_y_reverse(limits = c(0.0005, -0.015), breaks  = c(0, -0.015, -0.010, -0.005, 0)) 
  
ggsave(here("04-pdf", "k562_ms2_plci_allethers_20251002.pdf"), width = 4, height = 4)
ggsave(here("05-png", "k562_ms2_plcitrt_20250626b.png"), width = 4, height = 4)

# t-tests
plci_tot = plcidata2 %>% 
  #filter(time_gro == 6) %>%
  filter(scheme == "linreg") %>% 
  arrange(eid) %>% 
  # calc total PLCI
  group_by(eid, pres_gro, time_gro, scheme) %>% 
  summarise(plci = sum(plci)) %>% 
  ungroup()


# now summarise across replicates
plci_summary <- plci_tot %>%
  group_by(pres_gro, time_gro, scheme) %>%
  summarise(
    mean_plci = mean(plci, na.rm = TRUE),
    sd_plci   = sd(plci, na.rm = TRUE),
    sem_plci  = sd(plci, na.rm = TRUE) / sqrt(n()),
    n         = n(),
    .groups = "drop"
  )

# Control condition
control <- plci_summary %>%
  filter(pres_gro == 1, time_gro == 6) %>%
  select(mean_plci, sem_plci) %>%
  rename(control_mean = mean_plci, control_sem = sem_plci)

# Join control onto all groups
plci_diff <- plci_summary %>%
  mutate(tmp = 1) %>%
  left_join(control %>% mutate(tmp = 1), by = "tmp") %>%
  select(-tmp) %>%
  mutate(
    diff_mean = mean_plci - control_mean,
    diff_sem  = sqrt(sem_plci^2 + control_sem^2)  # error propagation
  )

# 6 h @ 1 vs. 250 bar
plci_tot %>% 
  filter(time_gro == 6) %>% 
  pivot_wider(names_from = "pres_gro", values_from = "plci") %>% 
  summarise(pairtest = t.test(
    `1`, `250`,
    alternative = "greater"
  ) %>% broom::tidy()) %>% 
  unnest(pairtest)

# 1 vs. 250 bar
# define a little Student-vs-Welch t-test wrapper
# checks equality of variances by F-test
# does NOT check normality of samples!
t.auto = function(x, y, ...){
  var_equal = (var.test(x, y) %>% broom::tidy() %>% .$p.value) > 0.05
  t.test(x, y, var.equal = var_equal, ...)
}

plci_tot %>% 
  filter(time_gro == 6) %>% 
  pivot_wider(names_from = "pres_gro", values_from = "plci") %>% 
  summarise(
    varitest = var.test(
      `1`, `250`,
    ) %>% broom::tidy(),
    pairtest = t.auto(
      `1`, `250`, 
      alternative = "greater"
    ) %>% broom::tidy()
  ) %>% 
  unnest(pairtest)
#####################################3
# t-tests
plci_tot = plcidata2 %>% 
  filter(scheme == "linreg") %>% 
  arrange(eid) %>% 
  # calc total PLCI
  group_by(eid, pres_gro, time_gro, scheme) %>% 
  summarise(plci = sum(plci)) %>% 
  ungroup()

# 6 h @ 1 vs. 250 bar
plci_tot %>% 
  #filter(time_gro == 6) %>% 
  pivot_wider(names_from = "pres_gro", values_from = "plci") %>% 
  summarise(pairtest = t.test(
    `1`, `250`,
    alternative = "greater"
  ) %>% broom::tidy()) %>% 
  unnest(pairtest)

# 1 vs. 250 bar
# define a little Student-vs-Welch t-test wrapper
# checks equality of variances by F-test
# does NOT check normality of samples!
t.auto = function(x, y, ...){
  var_equal = (var.test(x, y) %>% broom::tidy() %>% .$p.value) > 0.05
  t.test(x, y, var.equal = var_equal, ...)
}

plci_tot %>% 
  filter(time_gro == 6) %>% 
  pivot_wider(names_from = "pres_gro", values_from = "plci") %>% 
  summarise(
    #varitest = var.test(
    #  `1`, `250`,
    #) %>% broom::tidy(),
    pairtest = t.auto(
      `1`, `250`, 
      alternative = "greater"
    ) %>% broom::tidy()
  ) %>% 
  unnest(pairtest)



###########################################
# OK, try plotting the pooled data
plcidata2 %>% 
  filter(scheme == "linreg") %>% 
  #filter(time_gro == 6) %>%
  #filter(eid %in% c("K562 1b 6hrs 1", "K562 1b 6hrs 2", "K562 1b 6hrs 3", "K562 250b 18hrs 1", "K562 250b 18hrs 2", "K562 250b 18hrs 3")) %>%
  arrange(eid) %>% 
  group_by(eid, pres_gro, time_gro, class, scheme) %>% 
  summarise(plci = sum(plci)) %>% 
  # average by treatment and also get SEMs
  group_by(eid, pres_gro, time_gro, scheme) %>%
  mutate(plci_tot = ifelse(row_number() == 1, sum(plci), NA)) %>% 
  group_by(pres_gro, class, scheme) %>%
  summarize(
    plci = mean(plci),
    plci_sem = sd(plci_tot, na.rm = TRUE) / sqrt(n()),
    plci_tot = mean(plci_tot, na.rm = TRUE)
  ) %>% 
  arrange(pres_gro) %>% 
  mutate(
    trt = str_glue("{pres_gro} bar"),
    trt = trt %>% factor(levels = unique(.))
  ) %>% 
  gg_plcurv(
    darkmode = FALSE,
    #rows = scheme,
    cols = trt,
    y = plci 
  ) +
  # there are two errorbar colors
  # so you can make the downward-pointing side white-on-black if you like
  geom_errorbar(
    aes(
      x = 0.3, # hardcoded into gg_plcurv()
      ymin = plci_tot - plci_sem,
      ymax = plci_tot
    ),
    color = "black",
    #color = "white",
    width = 0
  ) +
  geom_errorbar(
    aes(
      x = 0.3, # hardcoded into gg_plcurv()
      ymin = plci_tot,
      ymax = plci_tot + plci_sem
    ),
    color = "black",
    width = 0
  ) +
  # make negative go upward!
  scale_y_reverse() +
  theme_pubr() +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    strip.text.x = element_text(angle = 90),
    legend.position = "right"
  ) +
  labs(
    #title = "K562 PL curvature by treatment",
    fill = "Lipid class",
    x = "Incubation pressure"
  )

ggsave(here("04-pdf", "k562_ms2_plcipress_sn2PUFAcorr_20250912.pdf"), width = 4, height = 4)
ggsave(here("05-png", "k562_ms2_plcipress_20250626a.png"), width = 4, height = 4)


#class composition
# OK, try plotting the pooled data
plcidata2 %>% 
  filter(scheme == "linreg") %>% 
  arrange(eid) %>% 
  group_by(eid, pres_gro, time_gro, class, scheme) %>% 
  summarise(frac_molar = sum(frac_molar)) %>% 
  # average by treatment and also get SEMs
  #group_by(eid, pres_gro, time_gro, scheme) %>%
  mutate(frac_tot = ifelse(row_number() == 1, sum(frac_molar), NA)) %>% 
  group_by(pres_gro, class, scheme) %>%
  summarize(
    frac_molar = mean(frac_molar),
    sem = sd(frac_tot, na.rm = TRUE) / sqrt(n()),
    tot = mean(frac_tot, na.rm = TRUE)
  ) %>% 
  arrange(pres_gro) %>% 
  mutate(
    trt = str_glue("{pres_gro} bar")#,
    #trt = trt %>% factor(levels = unique(.))
  ) %>% 
  gg_plcurv(
    darkmode = FALSE,
    #rows = scheme,
    cols = trt,
    y = frac_molar
  ) +
  # there are two errorbar colors
  # so you can make the downward-pointing side white-on-black if you like
  geom_errorbar(
    aes(
      x = 0.3, # hardcoded into gg_plcurv()
      ymin = frac_tot - sem,
      ymax = frac_tot
    ),
    color = "black",
    #color = "white",
    width = 0
  ) +
  geom_errorbar(
    aes(
      x = 0.3, # hardcoded into gg_plcurv()
      ymin = frac_tot,
      ymax = frac_tot + sem
    ),
    color = "black",
    width = 0
  ) +
  # make negative go upward!
  #scale_y_reverse() +
  theme_pubr() +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    strip.text.x = element_text(angle = 90),
    legend.position = "right"
  ) +
  labs(
    title = "K562 PL comp",
    fill = "Lipid class",
    x = "Incubation pressure"
  )


# Plot with the customized facet order
classes <- pldata %>% 
  arrange(eid) %>% 
  #filter(eid %in% c("K562 1b 6hrs 1","K562 1b 6hrs 2","K562 1b 6hrs 3","K562 250b 18hrs 1","K562 250b 18hrs 2", "K562 250b 18hrs 3" )) %>%
  group_by(pres_gro, time_gro, class, rep) %>% 
  summarise(frac_molar = sum(frac_molar), .groups = "drop_last") %>% 
  group_by(pres_gro, class) %>%
  summarise(
    mean_frac = mean(frac_molar, na.rm = TRUE),
    sem = sd(frac_molar, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
#%>%
ggplot() +
  #facet_wrap(~ pres_gro) +
  labs(x = "Pressure (bar)", y = "mole %") +
  scale_x_continuous(breaks = pldata %>% pull(pres_gro) %>% unique) +
  geom_col(
    data = classes, 
    aes(
      x = pres_gro,
      y = mean_frac*100,
      #fill = class
      fill = factor(class, levels = c("LPC", "LPE", "PS", "PC", "PI", "PE", "O-PC", "P-PC","O-PE", "P-PE" ))
    )
  ) +
  labs(fill = "Lipid") +
  theme_pubclean(base_size = 6) +
  scale_fill_manual(values = chroma_cl) +  
  #scale_fill_paletteer_d(("ggsci::category20_d3")) + #scale_fill_manual(values = chroma_cl) +
  #scale_fill_locuszoom() +
  theme(legend.position = "none") #+ 
theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())

ggsave(here("04-pdf", "k562_ms2_classcomp_allethers_20251002.pdf"), width = 2, height = 5.5, units = "cm")

#ethers#######################
##############
classes_subset <- plcidata2 %>%
  filter(class %in% c("P-PE", "O-PC", "P-PC", "O-PE")) %>%
  filter(pres_gro %in% c(1, 250)) 

# Plot
ggplot(classes_subset, aes(x = factor(pres_gro), y = mean_frac, fill = class)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_frac - sem, ymax = mean_frac + sem),
    position = position_dodge(width = 0.7),
    width = 0.3
  ) +
  labs(
    x = "Pressure (bar)",
    y = "Mole fraction",
    fill = "Lipid Class"
  ) +
  scale_x_discrete(labels = c("1", "250")) +
  theme_minimal(base_size = 10) +
  scale_fill_manual(values = chroma_cl) 
ggsave(here("04-pdf", "k562_PE_PC_ethers_20250626b.pdf"), width = 4, height = 4)

replicates <- pldata %>%
  #filter(time_gro == 6) %>%
  #filter(eid %in% c("K562 1b 6hrs 1","K562 1b 6hrs 2","K562 1b 6hrs 3","K562 250b 18hrs 1","K562 250b 18hrs 2", "K562 250b 18hrs 3" )) %>%
  group_by(eid, pres_gro, time_gro, class) %>%
  summarise(frac_molar = sum(frac_molar), .groups = "drop")
replicates_subset <- replicates %>%
  filter(class %in% c("P-PE", "O-PC", "P-PC", "O-PE")) %>%
  filter(pres_gro %in% c(1, 250))
library(dplyr)
library(broom)  # for tidy output
t_test_results <- replicates_subset %>%
  group_by(class) %>%
  summarise(
    t_test = list(t.test(frac_molar ~ pres_gro)),
    .groups = "drop"
  ) %>%
  mutate(tidy_results = map(t_test, broom::tidy)) %>%
  unnest(tidy_results)
t_test_results %>%
  select(class, estimate1, estimate2, statistic, p.value, conf.low, conf.high)


############
#combined ether vs ester
ether_sum <- pldata %>%
  filter(class %in% c("P-PE", "O-PC", "P-PC", "O-PE")) %>%
  #filter(time_gro == 6) %>%
  group_by(pres_gro, time_gro, eid) %>%
  summarise(ether_frac = sum(frac_molar), .groups = "drop")
ether_summary <- ether_sum %>%
  group_by(pres_gro) %>%
  summarise(
    mean_ether = mean(ether_frac),
    sem_ether = sd(ether_frac) / sqrt(n()),
    .groups = "drop"
  )
ggplot(ether_summary, aes(x = factor(pres_gro), y = mean_ether)) +
  geom_col(width = 0.6, fill = "steelblue") +
  geom_errorbar(
    aes(ymin = mean_ether - sem_ether, ymax = mean_ether + sem_ether),
    width = 0.3
  ) +
  labs(
    x = "Pressure (bar)",
    y = "Ether lipid molar fraction",
  ) +
  theme_minimal(base_size = 10)
ether_ttest <- t.test(ether_frac ~ pres_gro, data = ether_sum)
# View results
ether_ttest
ggsave(here("04-pdf", "k562_ethers_20250807.pdf"), width = 2, height = 4)



# Step 1: Get replicate-level molar fractions at 1 and 250 bar
lipid_reps <- pldata %>%
  filter(time_gro == 6) %>%
  filter(class %in% c("P-PE", "O-PC", "P-PC", "O-PE")) %>%
  filter(pres_gro %in% c(1, 250)) %>%
  group_by(class, pres_gro, rep) %>%
  summarise(frac_molar = sum(frac_molar), .groups = "drop") %>%
  pivot_wider(
    names_from = pres_gro,
    values_from = frac_molar,
    names_prefix = "bar_"
  ) %>%
  mutate(delta = bar_250 - bar_1)
lipid_delta_summary <- lipid_reps %>%
  group_by(class) %>%
  summarise(
    mean_delta = mean(delta, na.rm = TRUE),
    sem_delta = sd(delta, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

ggplot(lipid_delta_summary, aes(x = class, y = mean_delta, fill = class)) +
  geom_col() +
  geom_errorbar(
    aes(ymin = mean_delta - sem_delta, ymax = mean_delta + sem_delta),
    width = 0.3
  ) +
  labs(
        y = expression(Delta ~ "Molar fraction"),
  ) +
  theme_minimal(base_size = 12) +
  scale_fill_manual(values = chroma_cl) +  
  theme(legend.position = "none")





##############################################
#ether/ester ratio
#############################################
lipid_ratios <- pldata %>%
  filter(time_gro == 6) %>%
  filter(class %in% c("P-PE", "O-PC", "P-PC", "O-PE", "PE", "PC")) %>%
  filter(pres_gro %in% c(1, 250)) %>%
  group_by(pres_gro, time_gro, rep, class) %>%
  summarise(frac_molar = sum(frac_molar), .groups = "drop") %>%
  mutate(type = case_when(
    class %in% c("P-PE", "O-PC", "P-PC", "O-PE") ~ "ether",
    class %in% c("PE", "PC") ~ "ester"
  )) %>%
  group_by(pres_gro, time_gro, rep, type) %>%
  summarise(sum_frac = sum(frac_molar), .groups = "drop") %>%
  pivot_wider(names_from = type, values_from = sum_frac) %>%
  mutate(ratio = ether / ester)
ratio_summary <- lipid_ratios %>%
  group_by(pres_gro) %>%
  summarise(
    mean_ratio = mean(ratio),
    sem_ratio = sd(ratio) / sqrt(n()),
    .groups = "drop"
  )
ggplot(ratio_summary, aes(x = factor(pres_gro), y = mean_ratio)) +
  geom_col(width = 0.6, fill = "tomato") +
  geom_errorbar(
    aes(ymin = mean_ratio - sem_ratio, ymax = mean_ratio + sem_ratio),
    width = 0.3
  ) +
  labs(
    x = "Pressure (bar)",
    y = "Ether/Ester Lipid Ratio",
    title = "Ether/Ester Lipid Ratio at 1 and 250 bar"
  ) +
  theme_minimal(base_size = 12)
ether_ester_ttest <- t.test(ratio ~ pres_gro, data = lipid_ratios)
tidy(ether_ester_ttest)

# Filter for 6-hour growth time
ether_6h <- ether_sum %>%
  filter(time_gro == 6)

# Plot for 6-hour only
ggplot(ether_6h, aes(x = factor(pres_gro), y = ether_frac)) +
  geom_boxplot(
    fill = NA,
    color = "black",
    width = 0.6
  ) +
  labs(
    x = "Pressure (bar)",
    y = "Ether/Ester ratio"
  ) +
  scale_y_continuous(
    limits = c(0.20, 0.30),
    breaks = c(0.20, 0.25, 0.30)
  ) +
  theme_pubr(base_size = 6)
  ggsave(here("04-pdf", "k562_etheresterRatio_allethers_20251002.pdf"), width = 2.2, height = 2.4, units = "cm")


#total db comp

db_maj_summary <- plcidata %>% 
  filter(scheme == "linreg") %>%
  filter(time_gro == 6)  %>%
  #filter(total_DB >= 0 & total_DB <= 4) %>%
  arrange(eid) %>%
  group_by(eid, pres_gro, time_gro, dbonds, scheme) %>%
  summarize(db_frac = sum(frac_molar)) %>%
  group_by(pres_gro, time_gro, dbonds) %>% 
  summarize(
    avg_frac = mean(db_frac),
    sdev_frac = sd(db_frac),
    serr_frac = sdev_frac/sqrt(n()),
  )

mean_dbonds_pg <- db_maj_summary %>%
  filter(dbonds %in% 0:7) %>%
  group_by(pres_gro) %>%
  summarize(
    mean_db = weighted.mean(dbonds, avg_frac)
  )


db_maj_summary %>% 
  filter(dbonds %in% (0:7)) %>% ###double check filtering is with %in%
  ggplot() +
  #facet_wrap(~eid) + 
  geom_col(
    #data = pl_maj_summary, 
    position = position_dodge(),
    aes(
      x = dbonds,
      y = avg_frac*100,
      fill = as.factor(pres_gro)
    )
  ) + 
  #facet_grid(cols = vars(key_strain[strain])) + ###commented whis out
  labs(
    x = "Total Double Bonds",
    y = "mol % phospolipids", 
    fill = "Pressure (bar)"
  ) +
  theme_pubr(base_size = 8)+
  scale_fill_grey() +
  theme(legend.position = "right")+ 
  geom_errorbar(
    position = position_dodge(), 
    aes(
      x = dbonds,
      ymin = (avg_frac - serr_frac)*100,
      ymax = (avg_frac + serr_frac)*100,
      group = pres_gro
    )
  ) +
  theme(legend.position = 'none') 
ggsave(here("04-pdf", "k562_ms2_dbcomp.pdf"), width = 4, height = 4)

db_maj_summary %>% 
  filter(dbonds %in% 0:6) %>%
  ggplot() +
  geom_col(
    position = position_dodge(),
    aes(
      x = dbonds,
      y = avg_frac * 100,
      fill = as.factor(pres_gro)
    )
  ) + 
  geom_errorbar(
    position = position_dodge(), 
    aes(
      x = dbonds,
      ymin = (avg_frac - serr_frac) * 100,
      ymax = (avg_frac + serr_frac) * 100,
      group = pres_gro
    )
  ) +
  # Dashed vertical mean lines colored by pressure group
  geom_vline(
    data = mean_dbonds_pg,
    aes(
      xintercept = mean_db,
      color = as.factor(pres_gro),
      group = as.factor(pres_gro)
    ),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  labs(
    x = "Total Double Bonds",
    y = "mol % phospholipids", 
    fill = "Pressure (bar)",
    color = "Pressure (bar)"
  ) +
  theme_pubr(base_size = 8) +
  scale_fill_grey(start = 0.3, end = 0.7) +
  scale_color_grey(start = 0.3, end = 0.7) +
  theme(legend.position = "none")
ggsave(here("04-pdf", "k562_ms2_dbcomp_avg.pdf"), width = 4, height = 4)


#total length comp
len_maj_summary <- plcidata %>% 
  filter(scheme == "linreg") %>%
  filter(carbon > 27, carbon %% 2 == 0) %>%
  filter(time_gro == 6)  %>%
  #filter(total_DB >= 0 & total_DB <= 4) %>%
  arrange(eid) %>%
  group_by(eid, pres_gro, time_gro, carbon, scheme) %>%
  summarize(len_frac = sum(frac_molar)) %>%
  group_by(pres_gro, time_gro, carbon) %>% 
  summarize(
    avg_frac = mean(len_frac),
    sdev_frac = sd(len_frac),
    serr_frac = sdev_frac/sqrt(n()),
  )

mean_lengths_pg <- len_maj_summary %>%
  filter(carbon > 27, carbon %% 2 == 0) %>%
  group_by(pres_gro) %>%
  summarize(
    mean_len = weighted.mean(carbon, avg_frac)
  )


len_maj_summary %>% 
  #filter(time_gro == 6) %>%
  filter(carbon > 27) %>% ###double check filtering is with %in%
  filter(carbon %% 2 == 0) %>%
  ggplot() +
  #facet_wrap(~eid) + 
  geom_col(
    #data = pl_maj_summary, 
    position = position_dodge(),
    aes(
      x = carbon,
      y = avg_frac*100,
      fill = as.factor(pres_gro)
    )
  ) + 
  #facet_grid(cols = vars(key_strain[strain])) + ###commented whis out
  labs(
    x = "Total Length",
    y = "mol % phospolipids", 
    fill = "Pressure (bar)"
  ) +
  theme_pubr(base_size = 8)+
  scale_fill_grey() +
  theme(legend.position = "right")+ 
  geom_errorbar(
    position = position_dodge(), 
    aes(
      x = carbon,
      ymin = (avg_frac - serr_frac)*100,
      ymax = (avg_frac + serr_frac)*100,
      group = pres_gro
    )
  ) +
  theme(legend.position = 'none') +
  scale_x_continuous(limits = c(27, 45), breaks = seq(28, 46, 2))  # Adjusted to fit your OD values
ggsave(here("04-pdf", "k562_ms2_lengthcomp.pdf"), width = 4, height = 4)

len_maj_summary %>% 
  filter(carbon > 27, carbon %% 2 == 0) %>%
  ggplot() +
  geom_col(
    position = position_dodge(),
    aes(
      x = carbon,
      y = avg_frac * 100,
      fill = as.factor(pres_gro)
    )
  ) + 
  geom_errorbar(
    position = position_dodge(), 
    aes(
      x = carbon,
      ymin = (avg_frac - serr_frac) * 100,
      ymax = (avg_frac + serr_frac) * 100,
      group = pres_gro
    )
  ) +
  # Dashed vertical lines for mean total length
  geom_vline(
    data = mean_lengths_pg,
    aes(
      xintercept = mean_len,
      color = as.factor(pres_gro),
      group = as.factor(pres_gro)
    ),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  labs(
    x = "Total Length",
    y = "mol % phospholipids", 
    fill = "Pressure (bar)",
    color = "Pressure (bar)"
  ) +
  theme_pubr(base_size = 8) +
  scale_fill_grey(start = 0.3, end = 0.7) +
  scale_color_grey(start = 0.3, end = 0.7) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(27, 45), breaks = seq(28, 46, 2))
ggsave(here("04-pdf", "k562_ms2_lengthcomp_avg.pdf"), width = 4, height = 4)





















library(tidyverse)
library(ggpubr)
# lmapdata <- readr::read_tsv("lipidmaps_wmeta.tsv")   # optional if needed

# Define the three condition rows you want
# (rename columns here if your names differ)
# pres_gro: 1 or 250 ; time_gro: 6 or 18 ; eid: sample id ; class: lipid class
# carbon, dbonds should be total chain length and total DB across sn-1+sn-2.
# If you instead have per-chain columns (carbsn1, carbsn2, dbonsn1, dbonsn2) uncomment next mutate.
pldata <- lmapdata %>%
  mutate(
    group = case_when(
      pres_gro == 1   & time_gro == 6  ~ "1 bar, 6 h",
      pres_gro == 250 & time_gro == 6  ~ "250 bar, 6 h",
      pres_gro == 250 & time_gro == 18 ~ "250 bar, 18 h",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group)) %>%
  # keep only phospholipid headgroups you listed
  filter(str_detect(class, 'P') & !str_detect(class, 'Cer')) %>% 
  filter(class %in% c("LPC","LPE","PS","PC","O-PC","P-PC","PI","PE","P-PE","O-PE")) %>%
  group_by(eid) %>%
  mutate(frac_molar = frac_molar / sum(frac_molar)) %>%  # renormalize within each sample
  ungroup()
# If you have per-chain columns instead of totals, use:
# mutate(carbon = carbsn1 + carbsn2, dbonds = dbonsn1 + dbonsn2)

# Nice facet order
class_order <- c("PA","PC","PG","PE","PI","PS","LPC","LPE","O-PC","P-PC","O-PE","P-PE")
pldata <- pldata %>% mutate(class = factor(class, levels = intersect(class_order, unique(class))))
group_order <- c("1 bar, 6 h","250 bar, 6 h","250 bar, 18 h")
pldata <- pldata %>% mutate(group = factor(group, levels = group_order))

headgroup_df <- pldata %>%
  group_by(group, eid, class) %>%
  summarise(frac = sum(frac_molar), .groups = "drop") %>%     # class fraction per sample
  group_by(group, class) %>%
  summarise(frac_mean = mean(frac), .groups = "drop")

gg_headgroups <- headgroup_df %>%
  ggplot(aes(x = group, y = frac_mean*100, fill = class)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(scale = 1), expand = c(0,0)) +
  labs(x = NULL, y = "mole percent (%)", fill = "Headgroup") +
  theme_pubr() +
  theme(legend.position = "right")
gg_headgroups

unsat_df <- pldata %>%
  group_by(group, eid, class, total_DB = dbonds) %>%
  summarise(frac = sum(frac_molar), .groups = "drop") %>%
  group_by(group, eid, class) %>%
  mutate(frac_in_class = frac / sum(frac)) %>%       # normalize within class per sample
  ungroup() %>%
  group_by(group, class, total_DB) %>%
  summarise(frac_mean = mean(frac_in_class), .groups = "drop")

gg_unsat <- unsat_df %>%
  ggplot(aes(x = 1, y = frac_mean*100, fill = factor(total_DB))) +
  geom_col(width = 0.8) +
  facet_grid(rows = vars(group), cols = vars(class), switch = "y") +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0,100), expand = c(0,0)) +
  scale_x_continuous(breaks = NULL) +
  labs(x = "Pressure (bar)", y = "mole percent (%)", fill = "Total DB") +
  theme_pubr() +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey95", colour = NA),
    legend.position = "right"
  )
gg_unsat

length_df <- pldata %>%
  group_by(group, eid, class, total_len = carbon) %>%
  summarise(frac = sum(frac_molar), .groups = "drop") %>%
  group_by(group, eid, class) %>%
  mutate(frac_in_class = frac / sum(frac)) %>% 
  ungroup() %>%
  group_by(group, class, total_len) %>%
  summarise(frac_mean = mean(frac_in_class), .groups = "drop")

gg_length <- length_df %>%
  ggplot(aes(x = 1, y = frac_mean*100, fill = factor(total_len))) +
  geom_col(width = 0.8) +
  facet_grid(rows = vars(group), cols = vars(class), switch = "y") +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0,100), expand = c(0,0)) +
  scale_x_continuous(breaks = NULL) +
  labs(x = "Pressure (bar)", y = "mole percent (%)", fill = "Total length") +
  theme_pubr() +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey95", colour = NA),
    legend.position = "right"
  )
gg_length

library(tidyverse)
library(ggpubr)

baseline  <- "1 bar, 6 h"
compare   <- c("250 bar, 6 h", "250 bar, 18 h")

# Small helper to compute deltas after pivoting wide on `group`
make_deltas <- function(df, id_cols, value_col = "frac_mean") {
  wide <- df %>%
    pivot_wider(names_from = group, values_from = all_of(value_col))
  # Δ vs baseline (percentage points)
  for (cmp in compare) {
    wide[[paste0("delta_pp_", gsub("[ ,]", "_", cmp))]] <-
      100 * (wide[[cmp]] - wide[[baseline]])
  }
  wide
}


unsat_means <- pldata %>%
  filter(dbonds < 7) %>%
  group_by(group, eid, class, total_DB = dbonds) %>%
  summarise(frac = sum(frac_molar), .groups = "drop") %>%
  group_by(group, eid, class) %>%                      # normalize within class
  mutate(frac_in_class = frac / sum(frac)) %>%
  ungroup() %>%
  group_by(group, class, total_DB) %>%
  summarise(frac_mean = mean(frac_in_class), .groups = "drop")

unsat_deltas <- make_deltas(unsat_means, id_cols = c("class","total_DB"))

# Top movers per comparison (DB level inside each class)
top_unsat_250_6 <- unsat_deltas %>%
  arrange(desc(abs(delta_pp_250_bar__6_h))) %>%
  transmute(class, total_DB,
            `Δ pp (250 bar,6h – 1 bar,6h)` = delta_pp_250_bar__6_h) %>%
  slice_head(n = 15)

top_unsat_250_18 <- unsat_deltas %>%
  arrange(desc(abs(delta_pp_250_bar__18_h))) %>%
  transmute(class, total_DB,
            `Δ pp (250 bar,18h – 1 bar,6h)` = delta_pp_250_bar__18_h) %>%
  slice_head(n = 15)

# Heatmap per comparison
gg_unsat_delta <- unsat_deltas %>%
  select(class, total_DB, starts_with("delta_pp_")) %>%
  pivot_longer(starts_with("delta_pp_"), names_to = "comparison", values_to = "delta_pp") %>%
  mutate(comparison = recode(comparison,
                             delta_pp_250_bar__18_h = "250 bar, 18 h – 1 bar, 6 h"
  )) %>%
  ggplot(aes(x = factor(total_DB), y = class, fill = delta_pp)) +
  geom_tile() +
  facet_wrap(~ comparison) +
  scale_fill_gradient2(name = "Δ (pp)", limits = c(-20,20)) +
  labs(x = "Total DB", y = "Class", title = "Unsaturation shifts within class") +
  theme_pubr()
gg_unsat_delta

length_means <- pldata %>%
  filter(carbon > 30) %>%
  group_by(group, eid, class, total_len = carbon) %>%
  summarise(frac = sum(frac_molar), .groups = "drop") %>%
  group_by(group, eid, class) %>%
  mutate(frac_in_class = frac / sum(frac)) %>%
  ungroup() %>%
  group_by(group, class, total_len) %>%
  summarise(frac_mean = mean(frac_in_class), .groups = "drop")

length_deltas <- make_deltas(length_means, id_cols = c("class","total_len"))

top_length_250_6 <- length_deltas %>%
  arrange(desc(abs(delta_pp_250_bar__6_h))) %>%
  transmute(class, total_len,
            `Δ pp (250 bar,6h – 1 bar,6h)` = delta_pp_250_bar__6_h) %>%
  slice_head(n = 15)

top_length_250_18 <- length_deltas %>%
  arrange(desc(abs(delta_pp_250_bar__18_h))) %>%
  transmute(class, total_len,
            `Δ pp (250 bar,18h – 1 bar,6h)` = delta_pp_250_bar__18_h) %>%
  slice_head(n = 15)

gg_length_delta <- length_deltas %>%
  select(class, total_len, starts_with("delta_pp_")) %>%
  pivot_longer(starts_with("delta_pp_"), names_to = "comparison", values_to = "delta_pp") %>%
  mutate(comparison = recode(comparison,
                             delta_pp_250_bar__6_h  = "250 bar, 6 h – 1 bar, 6 h",
                             delta_pp_250_bar__18_h = "250 bar, 18 h – 1 bar, 6 h"
  )) %>%
  ggplot(aes(x = factor(total_len), y = class, fill = delta_pp)) +
  geom_tile() +
  facet_wrap(~ comparison) +
  scale_fill_gradient2(name = "Δ (pp)", midpoint = 0) +
  labs(x = "Total length", y = "Class", title = "Chain-length shifts within class") +
  theme_pubr()
gg_length_delta


library(tidyverse)
library(ggpubr)

## 1) Define the 3 experiment groups and clean up factors
pl_sn <- lmapdata %>%
  mutate(
    group = case_when(
      pres_gro == 1   & time_gro == 6  ~ "1 bar, 6 h",
      pres_gro == 250 & time_gro == 6  ~ "250 bar, 6 h",
      pres_gro == 250 & time_gro == 18 ~ "250 bar, 18 h",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group)) %>%
  # keep phospholipids you care about
  filter(str_detect(class, "P") & !str_detect(class, "Cer")) %>%
  filter(class %in% c("LPC","LPE","PS","PC","O-PC","P-PC","PI","PE","P-PE","O-PE")) %>%
  group_by(eid) %>%
  mutate(frac_molar = frac_molar / sum(frac_molar)) %>%  # renormalize per-sample to PL total
  ungroup()

# nice facet order (edit as you like)
class_order <- c("PA","PC","PG","PE","PI","PS","LPC","LPE","O-PC","P-PC","O-PE","P-PE")
pl_sn <- pl_sn %>%
  mutate(
    class = factor(class, levels = intersect(class_order, unique(class))),
    group = factor(group, levels = c("1 bar, 6 h","250 bar, 6 h","250 bar, 18 h"))
  )

## 2) Long format by sn position
sn1 <- pl_sn %>%
  transmute(eid, group, class, pos = "sn1",
            carbon = as.integer(round(carbsn1)),
            dbonds = as.integer(round(dbonsn1)),
            frac_molar)

sn2 <- pl_sn %>%
  transmute(eid, group, class, pos = "sn2",
            carbon = as.integer(round(carbsn2)),
            dbonds = as.integer(round(dbonsn2)),
            frac_molar)

sn_long <- bind_rows(sn1, sn2) %>%
  # drop chains that are missing (e.g., lysophospholipids will only contribute one sn)
  filter(!is.na(carbon), !is.na(dbonds))

## 3A) Within-class, within-position chain-length distributions
len_pos <- sn_long %>%
  group_by(group, eid, class, pos, carbon) %>%
  summarise(frac = sum(frac_molar), .groups = "drop") %>%
  group_by(group, eid, class, pos) %>%
  mutate(frac_in_class_pos = frac / sum(frac)) %>%   # normalize to 100% within class & position
  ungroup() %>%
  group_by(group, class, pos, carbon) %>%
  summarise(frac_mean = mean(frac_in_class_pos), .groups = "drop")

# order carbon bins if desired (e.g., 22–44)
len_levels <- sort(unique(len_pos$carbon))
len_pos <- len_pos %>% mutate(carbon = factor(carbon, levels = len_levels))

library(forcats)

p_len <- len_pos %>%
  filter(pos == "sn2",
         group %in% c("1 bar, 6 h", "250 bar, 18 h")) %>%
  # keep only even total lengths
  mutate(carbon_num = as.integer(as.character(carbon))) %>%
  filter(!is.na(carbon_num), carbon_num %% 2 == 0) %>%
  mutate(
    carbon = factor(carbon_num, levels = sort(unique(carbon_num))) # tidy legend order
  ) %>%
  ggplot(aes(x = pos, y = frac_mean*100, fill = carbon)) +
  geom_col(width = 0.8) +
  facet_grid(rows = vars(class), cols = vars(group), switch = "y") +
  scale_fill_carto_d(palette = "ag_Sunset", drop = TRUE) +  # drop odd levels from legend
  scale_y_continuous(limits = c(0,100), expand = c(0,0)) +
  labs(x = NULL, y = "mole percent within class (%)", fill = "sn2 length") +
  theme_pubr(base_size = 6) +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey95", colour = NA),
    legend.position = "right"
  )


## 3B) Within-class, within-position unsaturation distributions
unsat_pos <- sn_long %>%
  group_by(group, eid, class, pos, total_DB = dbonds) %>%
  summarise(frac = sum(frac_molar), .groups = "drop") %>%
  group_by(group, eid, class, pos) %>%
  mutate(frac_in_class_pos = frac / sum(frac)) %>%
  ungroup() %>%
  group_by(group, class, pos, total_DB) %>%
  summarise(frac_mean = mean(frac_in_class_pos), .groups = "drop")

unsat_levels <- sort(unique(unsat_pos$total_DB))
unsat_pos <- unsat_pos %>% mutate(total_DB = factor(total_DB, levels = unsat_levels))


library(rcartocolor)

p_unsat <- unsat_pos %>%
  filter(group %in% c("1 bar, 6 h", "250 bar, 18 h")) %>%
  filter(pos == "sn2") %>%
  ggplot(aes(x = pos, y = frac_mean*100, fill = total_DB)) +
  geom_col(width = 0.8) +
  #facet_grid(rows = vars(group), cols = vars(class), switch = "y") +
  facet_grid(rows = vars(class), cols = vars(group), switch = "y") +
  scale_fill_carto_d(palette = "ag_Sunset") +
  scale_y_continuous(limits = c(0,100), expand = c(0,0)) +
  labs(x = NULL, y = "mole percent within class (%)", fill = "sn2 DB") +
  theme_pubr(base_size = 6) +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey95", colour = NA),
    legend.position = "right"
  )

# Draw the plots
p_len
p_unsat
ggsave(here("04-pdf", "k562_sn2_length_stackbar.pdf"), width = 5, height = 18, units = "cm")




library(tidyverse)
library(ggpubr)

# -----------------------------
# 1) Prepare data & groups
# -----------------------------
pl_sn <- lmapdata %>%
  mutate(
    group = case_when(
      pres_gro == 1   & time_gro == 6 ~ "1 bar, 6 h",
      pres_gro == 250 & time_gro == 6 ~ "250 bar, 6 h",
      pres_gro == 250 & time_gro == 18 ~ "250 bar, 18 h"  #unused
    )
  ) %>%
  filter(group %in% c("1 bar, 6 h", "250 bar, 18 h")) %>%
  # keep phospholipids & renormalize to PL per sample (optional but recommended)
  filter(str_detect(class, "P") & !str_detect(class, "Cer")) %>%
  group_by(eid) %>%
  mutate(frac_molar = frac_molar / sum(frac_molar)) %>%
  ungroup()

# -----------------------------
# 2) Build sn2 tables (length + DB)
# -----------------------------
sn2 <- pl_sn %>%
  transmute(
    eid, group,
    len  = as.integer(round(carbsn2)),      # total carbons at sn2
    db   = as.integer(round(dbonsn2)),      # total double bonds at sn2
    w    = frac_molar                        # weight = mol fraction of species
  ) %>%
  filter(!is.na(len), !is.na(db))           # drop lysos that lack sn2

# (Optional) restrict to even lengths like your reference figure:
length_bins <- seq(0, 22, by = 2)
sn2 <- sn2 %>% filter(len %in% length_bins)

# -----------------------------
# 3) Distributions per replicate
# -----------------------------
# Length distribution (each bar = mean across replicates of the per-replicate % at that length)
len_rep <- sn2 %>%
  group_by(group, eid, len) %>%
  summarise(frac = sum(w), .groups = "drop") %>%
  group_by(group, eid) %>%
  mutate(frac = frac / sum(frac)) %>%                 # normalize within sn2 for that replicate
  ungroup()

# Fill missing length bins with 0 so SEM is defined
len_rep <- len_rep %>%
  group_by(group, eid) %>%
  tidyr::complete(len = length_bins, fill = list(frac = 0)) %>%
  ungroup()

# Unsaturation (DB) distribution
db_bins <- sort(unique(sn2$db))           # e.g., 0:6
db_rep <- sn2 %>%
  group_by(group, eid, db) %>%
  summarise(frac = sum(w), .groups = "drop") %>%
  group_by(group, eid) %>%
  mutate(frac = frac / sum(frac)) %>%
  ungroup() %>%
  group_by(group, eid) %>%
  tidyr::complete(db = db_bins, fill = list(frac = 0)) %>%
  ungroup()

# -----------------------------
# 4) Condition means, SEMs, and mean markers
# -----------------------------
sem <- function(x) sd(x)/sqrt(length(x))

len_sum <- len_rep %>%
  group_by(group, len) %>%
  summarise(mean = mean(frac), se = sem(frac), .groups = "drop")

db_sum <- db_rep %>%
  group_by(group, db) %>%
  summarise(mean = mean(frac), se = sem(frac), .groups = "drop")

# Weighted-mean markers (vertical lines)
len_mark <- sn2 %>%
  group_by(group, eid) %>%
  summarise(mu = sum(len * w) / sum(w), .groups = "drop") %>%
  group_by(group) %>%
  summarise(mu = mean(mu), .groups = "drop")

db_mark <- sn2 %>%
  group_by(group, eid) %>%
  summarise(mu = sum(db * w) / sum(w), .groups = "drop") %>%
  group_by(group) %>%
  summarise(mu = mean(mu), .groups = "drop")




# helper: Hedges' g (unbiased Cohen's d)
hedges_g <- function(x, y){
  n1 <- length(x); n2 <- length(y)
  s_p <- sqrt(((n1-1)*var(x) + (n2-1)*var(y)) / (n1+n2-2))
  d   <- (mean(x) - mean(y)) / s_p
  d * (1 - 3/(4*(n1+n2)-9))
}

# ============ sn2 LENGTH: per-bin tests ============
len_tests <- len_rep %>%
  group_by(len) %>%
  group_modify(~{
    x <- .x$frac[.x$group == "1 bar, 6 h"]
    y <- .x$frac[.x$group == "250 bar, 6 h"]
    tt <- t.test(y, x, var.equal = FALSE)  # (250,18) − (1,6)
    tibble(
      n_1bar = length(x),
      n_250  = length(y),
      mean_1bar = mean(x)*100,
      mean_250  = mean(y)*100,
      diff_250_minus_1 = (mean(y) - mean(x))*100,    # percentage points
      p_t   = tt$p.value,
      p_wcx = tryCatch(wilcox.test(y, x)$p.value, error = function(e) NA_real_),
      hedges_g = hedges_g(y, x)
    )
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_t, method = "BH"),
         sig = case_when(
           p_adj < 0.001 ~ "***",
           p_adj < 0.01  ~ "**",
           p_adj < 0.05  ~ "*",
           TRUE ~ ""
         )) %>%
  arrange(p_adj)

len_tests


# ============ sn2 UNSATURATION (DB): per-bin tests ============
db_tests <- db_rep %>%
  group_by(db) %>%
  group_modify(~{
    x <- .x$frac[.x$group == "1 bar, 6 h"]
    y <- .x$frac[.x$group == "250 bar, 6 h"]
    tt <- t.test(y, x, var.equal = FALSE)
    tibble(
      n_1bar = length(x),
      n_250  = length(y),
      mean_1bar = mean(x)*100,
      mean_250  = mean(y)*100,
      diff_250_minus_1 = (mean(y) - mean(x))*100,
      p_t   = tt$p.value,
      p_wcx = tryCatch(wilcox.test(y, x)$p.value, error = function(e) NA_real_),
      hedges_g = hedges_g(y, x)
    )
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_t, method = "BH"),
         sig = case_when(
           p_adj < 0.001 ~ "***",
           p_adj < 0.01  ~ "**",
           p_adj < 0.05  ~ "*",
           TRUE ~ ""
         )) %>%
  arrange(p_adj)

db_tests



# --------------------------------------------
# 1) Build sn1 + sn2 long table
# --------------------------------------------
sn_long <- bind_rows(
  pl_sn %>%
    transmute(eid, group, pos = "sn1",
              len = as.integer(round(carbsn1)),
              db  = as.integer(round(dbonsn1)),
              w   = frac_molar),
  pl_sn %>%
    transmute(eid, group, pos = "sn2",
              len = as.integer(round(carbsn2)),
              db  = as.integer(round(dbonsn2)),
              w   = frac_molar)
) %>%
  filter(!is.na(len), !is.na(db))

# keep only even lengths? (set to FALSE to keep all)
even_only <- TRUE
if (even_only) sn_long <- sn_long %>% filter(len %% 2 == 0)

# define contiguous bin sets so missing bins can be filled with 0s
len_levels <- {
  mn <- min(sn_long$len); mx <- max(sn_long$len)
  if (even_only) seq(ifelse(mn %% 2 == 0, mn, mn+1),
                     ifelse(mx %% 2 == 0, mx, mx-1), by = 2)
  else sort(unique(sn_long$len))
}
db_levels  <- sort(unique(sn_long$db))

# --------------------------------------------
# 2) Per-replicate distributions (normalize within pos)
# --------------------------------------------
# Length
len_rep_pos <- sn_long %>%
  group_by(group, eid, pos, len) %>%
  summarise(frac = sum(w), .groups = "drop") %>%
  group_by(group, eid, pos) %>%
  mutate(frac = frac / sum(frac)) %>%         # distribution within sn-position
  ungroup() %>%
  group_by(group, eid, pos) %>%               # fill missing bins with 0
  complete(len = len_levels, fill = list(frac = 0)) %>%
  ungroup()

# Double bonds
db_rep_pos <- sn_long %>%
  group_by(group, eid, pos, db) %>%
  summarise(frac = sum(w), .groups = "drop") %>%
  group_by(group, eid, pos) %>%
  mutate(frac = frac / sum(frac)) %>%
  ungroup() %>%
  group_by(group, eid, pos) %>%
  complete(db = db_levels, fill = list(frac = 0)) %>%
  ungroup()

# --------------------------------------------
# 3) Per-bin tests (Welch t, Wilcoxon, BH FDR) + effect size
# --------------------------------------------
hedges_g <- function(x, y){
  n1 <- length(x); n2 <- length(y)
  sp <- sqrt(((n1-1)*var(x) + (n2-1)*var(y)) / (n1+n2-2))
  d  <- (mean(y) - mean(x)) / sp                # (250bar,18h) – (1bar,6h)
  d * (1 - 3/(4*(n1+n2) - 9))
}

# Length tests, by sn-position
len_tests_pos <- len_rep_pos %>%
  group_by(pos, len) %>%
  group_modify(~{
    x <- .x$frac[.x$group == "1 bar, 6 h"]
    y <- .x$frac[.x$group == "250 bar, 18 h"]
    tt <- t.test(y, x, var.equal = FALSE)
    tibble(
      n_1bar = length(x), n_250 = length(y),
      mean_1bar = mean(x)*100, mean_250 = mean(y)*100,
      diff_250_minus_1 = (mean(y) - mean(x))*100,
      p_t = tt$p.value,
      p_wilcox = tryCatch(wilcox.test(y, x)$p.value, error = function(e) NA_real_),
      hedges_g = hedges_g(x, y)
    )
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_t, method = "BH"),
         sig = case_when(p_adj < 0.001 ~ "***",
                         p_adj < 0.01  ~ "**",
                         p_adj < 0.05  ~ "*",
                         TRUE ~ "")) %>%
  arrange(pos, p_adj)

# DB tests, by sn-position
db_tests_pos <- db_rep_pos %>%
  group_by(pos, db) %>%
  group_modify(~{
    x <- .x$frac[.x$group == "1 bar, 6 h"]
    y <- .x$frac[.x$group == "250 bar, 18 h"]
    tt <- t.test(y, x, var.equal = FALSE)
    tibble(
      n_1bar = length(x), n_250 = length(y),
      mean_1bar = mean(x)*100, mean_250 = mean(y)*100,
      diff_250_minus_1 = (mean(y) - mean(x))*100,
      p_t = tt$p.value,
      p_wilcox = tryCatch(wilcox.test(y, x)$p.value, error = function(e) NA_real_),
      hedges_g = hedges_g(x, y)
    )
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_t, method = "BH"),
         sig = case_when(p_adj < 0.001 ~ "***",
                         p_adj < 0.01  ~ "**",
                         p_adj < 0.05  ~ "*",
                         TRUE ~ "")) %>%
  arrange(pos, p_adj)

len_tests_pos
db_tests_pos


# Greys for bars/lines
custom_greys <- c("1 bar, 6 h" = "#222222", "250 bar, 18 h" = "#A9A9A9")

# -----------------------------
# 5A) Plot: sn2 total length distribution
# -----------------------------
pd <- position_dodge(width = 0.9)

p_len <- ggplot(len_sum, aes(x = len, y = mean*100, fill = group)) +
  geom_col(position = pd, width = 0.85) +
  geom_errorbar(
    aes(ymin = (mean - se)*100, ymax = (mean + se)*100),
    position = pd, width = 0.85, linewidth = 0.3
  ) +
  # mean markers (one dashed line per condition)
  geom_vline(data = len_mark, aes(xintercept = mu, color = group),
             linetype = "dashed", linewidth = 0.6, alpha = 0.9, show.legend = FALSE) +
  scale_fill_manual(values = custom_greys) +
  scale_color_manual(values = custom_greys) +
  scale_x_continuous(limits = c(13, 23), breaks = c(14, 16, 18, 20, 22, 24)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Total Length (sn-2)", y = "mol % phospholipids", fill = NULL) +
  theme_pubr(base_size = 6) +
  theme(legend.position = "none")

# -----------------------------
# 5B) Plot: sn2 total double bonds distribution
# -----------------------------
p_db <- ggplot(db_sum, aes(x = db, y = mean*100, fill = group)) +
  geom_col(position = pd, width = 0.85) +
  geom_errorbar(
    aes(ymin = (mean - se)*100, ymax = (mean + se)*100),
    position = pd, width = 0.85, linewidth = 0.3) +
  geom_vline(data = db_mark, aes(xintercept = mu, color = group),
             linetype = "dashed", linewidth = 0.6, alpha = 0.9, show.legend = FALSE) +
  scale_fill_manual(values = custom_greys) +
  scale_color_manual(values = custom_greys) +
  scale_x_continuous(breaks = db_bins) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Total Double Bonds (sn-2)", y = "mol % phospholipids", fill = NULL) +
  theme_pubr(base_size = 6) +
  theme(legend.position = "none")

p_len
ggsave(here("04-pdf", "k562_sn2_length.pdf"), width = 3, height = 3, units = "cm")

p_db
ggsave(here("04-pdf", "k562_sn2_db.pdf"), width = 3, height = 3, units = "cm")


#sn2 total mean significance tests
library(tidyverse)
library(broom)

# 1) Build the per–sn2 table (make sure these column names match your data)
pl_sn <- lmapdata %>%
  mutate(
    group = case_when(
      pres_gro == 1   & time_gro == 6 ~ "1 bar, 6 h",
      pres_gro == 250 & time_gro == 6 ~ "250 bar, 6 h",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group)) %>%
  filter(str_detect(class, "P") & !str_detect(class, "Cer")) %>%
  group_by(eid) %>% mutate(frac_molar = frac_molar / sum(frac_molar)) %>% ungroup()

sn2 <- pl_sn %>%
  transmute(
    eid, group, class,
    len = as.integer(round(carbsn2)),       # total carbons at sn-2
    db  = as.integer(round(dbonsn2)),       # total double bonds at sn-2
    w   = frac_molar                        # weight = mol fraction
  ) %>%
  filter(!is.na(len), !is.na(db))

# 2) Replicate-level weighted means (this creates mu_len and mu_db)
rep_means <- sn2 %>%
  group_by(group, eid) %>%
  summarise(
    mu_len = weighted.mean(len, w),
    mu_db  = weighted.mean(db,  w),
    .groups = "drop"
  )

# quick sanity check (optional): head(rep_means)

# 3) Significance tests (Welch t-test; add Wilcoxon if you like)
len_t <- t.test(mu_len ~ group, data = rep_means, var.equal = FALSE)
db_t  <- t.test(mu_db  ~ group, data = rep_means, var.equal = FALSE)

tidy(len_t)
tidy(db_t)

by_class <- sn2 %>%
  group_by(class, group, eid) %>%
  summarise(
    mu_len = weighted.mean(len, w),
    mu_db  = weighted.mean(db,  w),
    .groups = "drop"
  )

class_tests <- by_class %>%
  group_by(class) %>%
  group_modify(~{
    tibble(
      class = unique(.x$class),
      len_p = if (n_distinct(.x$group) == 2) t.test(mu_len ~ group, data = .x)$p.value else NA_real_,
      db_p  = if (n_distinct(.x$group) == 2) t.test(mu_db  ~ group, data = .x)$p.value else NA_real_
    )
  }) %>% ungroup() %>%
  mutate(p_adj_len = p.adjust(len_p, method = "BH"),
         p_adj_db  = p.adjust(db_p,  method = "BH"))

class_tests


library(tidyverse)
library(ggpubr)

# ---------- 1) Prepare data ----------
even_only <- TRUE                       # set FALSE to include odd lengths too

pl_sn <- lmapdata %>%
  mutate(
    group = case_when(
      pres_gro == 1   & time_gro == 6  ~ "1 bar, 6 h",
      pres_gro == 250 & time_gro == 6  ~ "250 bar, 6 h",
      pres_gro == 250 & time_gro == 18 ~ "250 bar, 18 h",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group)) %>%
  filter(str_detect(class, "P") & !str_detect(class, "PG")) %>%   # phospholipids only
  group_by(eid) %>% mutate(frac_molar = frac_molar / sum(frac_molar)) %>% ungroup()

sn2 <- pl_sn %>%
  transmute(
    eid, group, class,
    len = as.integer(round(carbsn2)),            # total carbons at sn-2
    db  = as.integer(round(dbonsn2)),            # total double bonds at sn-2
    w   = frac_molar
  ) %>%
  filter(!is.na(len), !is.na(db))

# bin sets
len_bins <- if (even_only) seq(0, 22, by = 2) else sort(unique(sn2$len))
db_bins  <- sort(unique(sn2$db))

# ---------- 2) Per-replicate sn2 distributions, normalized within class ----------
# Length
len_rep <- sn2 %>%
  group_by(group, eid, class, len) %>%
  summarise(frac = sum(w), .groups = "drop") %>%
  group_by(group, eid, class) %>%
  mutate(frac = frac / sum(frac)) %>%                # within-class normalization
  ungroup() %>%
  group_by(group, eid, class) %>%                    # fill missing bins
  complete(len = len_bins, fill = list(frac = 0)) %>%
  ungroup()

# Unsaturation (DB)
db_rep <- sn2 %>%
  group_by(group, eid, class, db) %>%
  summarise(frac = sum(w), .groups = "drop") %>%
  group_by(group, eid, class) %>%
  mutate(frac = frac / sum(frac)) %>%
  ungroup() %>%
  group_by(group, eid, class) %>%
  complete(db = db_bins, fill = list(frac = 0)) %>%
  ungroup()

# ---------- 3) Condition means (percent within class) ----------
len_mean <- len_rep %>%
  group_by(group, class, len) %>%
  summarise(pct = mean(frac)*100, .groups = "drop")

db_mean <- db_rep %>%
  group_by(group, class, db) %>%
  summarise(pct = mean(frac)*100, .groups = "drop")

# Length deltas
len_delta <- len_mean %>%
  filter(len > 12) %>%
  pivot_wider(names_from = group, values_from = pct) %>%
  mutate(`Δ 250-1 (pp)` = `250 bar, 18 h` - `1 bar, 6 h`) %>%
  select(class, len, `Δ 250-1 (pp)`)

p_len_delta <- ggplot(len_delta, aes(x = factor(len), y = class, fill = `Δ 250-1 (pp)`)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0, name = "Δ (pp)") +
  labs(x = "sn-2 total length", y = NULL, title = "Change: 250 bar, 18 h – 1 bar, 6 h") +
  theme_pubr()

# Unsaturation deltas
db_delta <- db_mean %>%
  pivot_wider(names_from = group, values_from = pct) %>%
  mutate(`Δ 250-1 (pp)` = `250 bar, 18 h` - `1 bar, 6 h`) %>%
  select(class, db, `Δ 250-1 (pp)`)

p_db_delta <- ggplot(db_delta, aes(x = factor(db), y = class, fill = `Δ 250-1 (pp)`)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "darkblue", mid = "white", high = "darkred",
    midpoint = 0, limits = c(-15, 15), name = "Δ mol %") +
  labs(x = "sn-2 double bonds", y = NULL)+
  theme_pubr(base_size = 6) +
  theme(legend.position = "none")

p_db_delta
ggsave(here("04-pdf", "k562_heatmap_allether_sn2_db.pdf"), width = 3, height = 3, units = "cm")


# starting from your len_delta data frame
# (columns: class, len, `Δ 250-1 (pp)` or similar)

# keep only even sn-2 lengths and make them the only x levels
len_delta_even <- len_delta %>% 
  filter(len > 14, len < 24) %>%
  mutate(len = as.integer(as.character(len))) %>%          # ensure numeric
  filter(!is.na(len), len %% 2 == 0)

# define contiguous even levels to avoid “skips”
min_even <- ifelse(min(len_delta_even$len) %% 2 == 0,
                   min(len_delta_even$len), min(len_delta_even$len) + 1)
max_even <- ifelse(max(len_delta_even$len) %% 2 == 0,
                   max(len_delta_even$len), max(len_delta_even$len) - 1)
even_levels <- seq(min_even, max_even, by = 2)

len_delta_even <- len_delta_even %>%
  mutate(len = factor(len, levels = even_levels, ordered = TRUE))

library(scales)  # for squish

p_len_delta <- ggplot(len_delta_even,
                      aes(x = len, y = class, fill = `Δ 250-1 (pp)`)) +
  geom_tile() +
  scale_x_discrete(drop = FALSE) +
  scale_fill_gradient2(
    low = "darkblue", mid = "white", high = "darkred",
    midpoint = 0, limits = c(-20, 20), oob = squish,
    name = "\u0394 mol %"           # legend title, no "pp"
  ) +
  labs(x = "sn-2 length", y = NULL) +
  theme_pubr(base_size = 6) +
  theme(legend.position="none")


p_len_delta
ggsave(here("04-pdf", "k562_heatmap_allether_length.pdf"), width = 3, height = 3, units = "cm")


library(tidyverse)
library(ggpubr)
library(broom)

#---------------------------------------------
# Choose the two groups to compare
#   (switch to "250 bar, 18 h" if you prefer)
groups_keep <- c("1 bar, 6 h", "250 bar, 18 h")

# Base PL table -> keep only those groups and renormalize per sample
pl_sn1 <- lmapdata %>%
  mutate(group = case_when(
    pres_gro == 1   & time_gro == 6  ~ "1 bar, 6 h",
    pres_gro == 250 & time_gro == 6  ~ "250 bar, 6 h",
    pres_gro == 250 & time_gro == 18 ~ "250 bar, 18 h",
    TRUE ~ NA_character_
  )) %>%
  filter(group %in% groups_keep) %>%
  filter(str_detect(class, "P") & !str_detect(class, "Cer")) %>%
  group_by(eid) %>% mutate(frac_molar = frac_molar / sum(frac_molar)) %>% ungroup()

# sn-1 table
sn1 <- pl_sn1 %>%
  transmute(eid, group,
            len = as.integer(round(carbsn1)),
            db  = as.integer(round(dbonsn1)),
            w   = frac_molar) %>%
  filter(!is.na(len), !is.na(db))

# Keep only even sn-1 lengths (like your example)
even_only <- TRUE
if (even_only) sn1 <- sn1 %>% filter(len > 12) %>% filter(len < 24) %>%filter(len %% 2 == 0) 

# Define contiguous bins
len_bins <- {
  mn <- min(sn1$len); mx <- max(sn1$len)
  if (even_only) seq(ifelse(mn %% 2 == 0, mn, mn+1),
                     ifelse(mx %% 2 == 0, mx, mx-1), by = 2)
  else sort(unique(sn1$len))
}
db_bins <- sort(unique(sn1$db))

#---------------------------------------------
# Distributions per replicate (normalize within sn-1)
len_rep1 <- sn1 %>%
  group_by(group, eid, len) %>% summarise(frac = sum(w), .groups = "drop") %>%
  group_by(group, eid) %>% mutate(frac = frac / sum(frac)) %>% ungroup() %>%
  group_by(group, eid) %>% complete(len = len_bins, fill = list(frac = 0)) %>% ungroup()

db_rep1 <- sn1 %>%
  group_by(group, eid, db) %>% summarise(frac = sum(w), .groups = "drop") %>%
  group_by(group, eid) %>% mutate(frac = frac / sum(frac)) %>% ungroup() %>%
  group_by(group, eid) %>% complete(db = db_bins, fill = list(frac = 0)) %>% ungroup()

sem <- function(x) sd(x)/sqrt(length(x))

len_sum1 <- len_rep1 %>%
  group_by(group, len) %>%
  summarise(mean = mean(frac), se = sem(frac), .groups = "drop")

db_sum1 <- db_rep1 %>%
  group_by(group, db) %>%
  summarise(mean = mean(frac), se = sem(frac), .groups = "drop")

# Mean markers (replicate-weighted means, then averaged across replicates)
len_mark1 <- sn1 %>%
  group_by(group, eid) %>% summarise(mu = sum(len*w)/sum(w), .groups = "drop") %>%
  group_by(group) %>% summarise(mu = mean(mu), .groups = "drop")

db_mark1 <- sn1 %>%
  group_by(group, eid) %>% summarise(mu = sum(db*w)/sum(w), .groups = "drop") %>%
  group_by(group) %>% summarise(mu = mean(mu), .groups = "drop")

#---------------------------------------------
# Per-bin significance (Welch t, BH FDR) for stars
hedges_g <- function(x, y){
  n1 <- length(x); n2 <- length(y)
  sp <- sqrt(((n1-1)*var(x) + (n2-1)*var(y)) / (n1+n2-2))
  d  <- (mean(y) - mean(x)) / sp
  d * (1 - 3/(4*(n1+n2) - 9))
}

len_tests1 <- len_rep1 %>%
  group_by(len) %>%
  group_modify(~{
    x <- .x$frac[.x$group == groups_keep[1]]
    y <- .x$frac[.x$group == groups_keep[2]]
    tt <- t.test(y, x, var.equal = FALSE)
    tibble(p = tt$p.value,
           diff_pp = (mean(y) - mean(x))*100,
           g = hedges_g(x, y))
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p, method = "BH"),
         sig = case_when(p_adj < 0.001 ~ "***",
                         p_adj < 0.01  ~ "**",
                         p_adj < 0.05  ~ "*",
                         TRUE ~ ""))

db_tests1 <- db_rep1 %>%
  group_by(db) %>%
  group_modify(~{
    x <- .x$frac[.x$group == groups_keep[1]]
    y <- .x$frac[.x$group == groups_keep[2]]
    tt <- t.test(y, x, var.equal = FALSE)
    tibble(p = tt$p.value,
           diff_pp = (mean(y) - mean(x))*100,
           g = hedges_g(x, y))
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p, method = "BH"),
         sig = case_when(p_adj < 0.001 ~ "***",
                         p_adj < 0.01  ~ "**",
                         p_adj < 0.05  ~ "*",
                         TRUE ~ ""))

# y-positions for stars (above tallest bar at each bin)
len_ann1 <- len_sum1 %>%
  group_by(len) %>% summarise(y = max((mean + se)*100) + 1, .groups = "drop") %>%
  left_join(len_tests1 %>% select(len, sig), by = "len")

db_ann1 <- db_sum1 %>%
  group_by(db) %>% summarise(y = max((mean + se)*100) + 1, .groups = "drop") %>%
  left_join(db_tests1 %>% select(db, sig), by = "db")

#---------------------------------------------
# Plotting (style matches your sn-2 figure)
custom_greys <- c("1 bar, 6 h" = "#555555", "250 bar, 18 h" = "#A9A9A9")
pd <- position_dodge(width = 0.9)

# A) sn-1 total length distribution
p_sn1_len <- ggplot(len_sum1, aes(x = len, y = mean*100, fill = group)) +
  geom_col(position = pd, width = 0.85) +
  geom_errorbar(aes(ymin = (mean - se)*100, ymax = (mean + se)*100),
                position = pd, width = 0.85, linewidth = 0.3) +
  geom_vline(data = len_mark1, aes(xintercept = mu, color = group),
             linetype = "dashed", linewidth = 0.6, alpha = 0.9, show.legend = FALSE) +
  #geom_text(data = len_ann1, aes(x = len, y = y, label = sig), size = 3, vjust = 0) +
  scale_fill_manual(values = custom_greys, breaks = groups_keep, name = NULL) +
  scale_color_manual(values = custom_greys, guide = "none") +
  scale_x_continuous(limits = c(13, 23), breaks = len_bins) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "sn-1 length", y = "mol % phospholipids") +
  theme_pubr(base_size = 5) +
  theme(legend.position = "none")


# B) sn-1 total double bonds distribution
p_sn1_db <- ggplot(db_sum1, aes(x = db, y = mean*100, fill = group)) +
  geom_col(position = pd, width = 0.85) +
  geom_errorbar(aes(ymin = (mean - se)*100, ymax = (mean + se)*100),
                position = pd, width = 0.85, linewidth = 0.3) +
  geom_vline(data = db_mark1, aes(xintercept = mu, color = group),
             linetype = "dashed", linewidth = 0.6, alpha = 0.9, show.legend = FALSE) +
  #geom_text(data = db_ann1, aes(x = db, y = y, label = sig), size = 3, vjust = 0) +
  scale_fill_manual(values = custom_greys, breaks = groups_keep, name = NULL) +
  scale_color_manual(values = custom_greys, guide = "none") +
  scale_x_continuous(breaks = db_bins) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "sn-1 double bonds", y = "mol % phospholipids") +
  theme_pubr(base_size = 5) +
  theme(legend.position = "none")

p_sn1_len
ggsave(here("04-pdf", "k562_sn1_length.pdf"), width = 3, height = 3, units = "cm")
p_sn1_db
ggsave(here("04-pdf", "k562_sn1_db.pdf"), width = 3, height = 3, units = "cm")




library(tidyverse)
library(ggpubr)
library(broom)

# --------------------------------------------
# Choose the two conditions to compare
# --------------------------------------------
groups_keep <- c("1 bar, 6 h", "250 bar, 18 h")   # or c("1 bar, 6 h","250 bar, 18 h")

# --------------------------------------------
# Build PL table: keep groups, PL only, renormalize per sample
# --------------------------------------------
pl_tot <- lmapdata %>%
  mutate(group = case_when(
    pres_gro == 1   & time_gro == 6  ~ "1 bar, 6 h",
    pres_gro == 250 & time_gro == 6  ~ "250 bar, 6 h",
    pres_gro == 250 & time_gro == 18 ~ "250 bar, 18 h",
    TRUE ~ NA_character_
  )) %>%
  filter(group %in% groups_keep) %>%
  filter(str_detect(class, "P") & !str_detect(class, "Cer")) %>%
  group_by(eid) %>% mutate(frac_molar = frac_molar / sum(frac_molar)) %>% ungroup()

# --------------------------------------------
# Total-length / total-DB table (sn1+sn2 combined)
# --------------------------------------------
tot <- pl_tot %>%
  transmute(
    eid, group,
    len = as.integer(round(carbon)),   # TOTAL carbons (sn1 + sn2)
    db  = as.integer(round(dbonds)),   # TOTAL double bonds (sn1 + sn2)
    w   = frac_molar                   # weight = mol fraction
  ) %>%
  filter(!is.na(len), !is.na(db))

# Keep only even total lengths (set to FALSE to include odds)
even_only <- TRUE
if (even_only) tot <- tot %>% filter(len %% 2 == 0)

# Define contiguous bins so there are no gaps on the x-axis
len_bins <- {
  mn <- min(tot$len); mx <- max(tot$len)
  if (even_only) seq(ifelse(mn %% 2 == 0, mn, mn+1),
                     ifelse(mx %% 2 == 0, mx, mx-1), by = 2)
  else sort(unique(tot$len))
}
db_bins <- sort(unique(tot$db))

# --------------------------------------------
# Per-replicate distributions (normalize within TOTAL set)
# --------------------------------------------
len_rep <- tot %>%
  group_by(group, eid, len) %>% summarise(frac = sum(w), .groups = "drop") %>%
  group_by(group, eid) %>% mutate(frac = frac / sum(frac)) %>% ungroup() %>%
  group_by(group, eid) %>% complete(len = len_bins, fill = list(frac = 0)) %>% ungroup()

db_rep <- tot %>%
  group_by(group, eid, db) %>% summarise(frac = sum(w), .groups = "drop") %>%
  group_by(group, eid) %>% mutate(frac = frac / sum(frac)) %>% ungroup() %>%
  group_by(group, eid) %>% complete(db = db_bins, fill = list(frac = 0)) %>% ungroup()

sem <- function(x) sd(x)/sqrt(length(x))

len_sum <- len_rep %>% group_by(group, len) %>% summarise(mean = mean(frac), se = sem(frac), .groups = "drop")
db_sum  <- db_rep  %>% group_by(group, db ) %>% summarise(mean = mean(frac), se = sem(frac), .groups = "drop")

# Replicate-weighted mean markers (dashed lines)
len_mark <- tot %>% group_by(group, eid) %>% summarise(mu = sum(len*w)/sum(w), .groups = "drop") %>%
  group_by(group) %>% summarise(mu = mean(mu), .groups = "drop")
db_mark  <- tot %>% group_by(group, eid) %>% summarise(mu = sum(db*w)/sum(w),  .groups = "drop") %>%
  group_by(group) %>% summarise(mu = mean(mu), .groups = "drop")

# --------------------------------------------
# Per-bin significance with BH FDR (for stars)
# --------------------------------------------
hedges_g <- function(x, y){
  n1 <- length(x); n2 <- length(y)
  sp <- sqrt(((n1-1)*var(x) + (n2-1)*var(y)) / (n1+n2-2))
  d  <- (mean(y) - mean(x)) / sp
  d * (1 - 3/(4*(n1+n2) - 9))
}

len_tests <- len_rep %>%
  group_by(len) %>%
  group_modify(~{
    x <- .x$frac[.x$group == groups_keep[1]]
    y <- .x$frac[.x$group == groups_keep[2]]
    tt <- t.test(y, x, var.equal = FALSE)
    tibble(p = tt$p.value, diff_pp = (mean(y) - mean(x))*100, g = hedges_g(x, y))
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p, method = "BH"),
         sig = case_when(p_adj < 0.001 ~ "***",
                         p_adj < 0.01  ~ "**",
                         p_adj < 0.05  ~ "*",
                         TRUE ~ ""))

db_tests <- db_rep %>%
  group_by(db) %>%
  group_modify(~{
    x <- .x$frac[.x$group == groups_keep[1]]
    y <- .x$frac[.x$group == groups_keep[2]]
    tt <- t.test(y, x, var.equal = FALSE)
    tibble(p = tt$p.value, diff_pp = (mean(y) - mean(x))*100, g = hedges_g(x, y))
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p, method = "BH"),
         sig = case_when(p_adj < 0.001 ~ "***",
                         p_adj < 0.01  ~ "**",
                         p_adj < 0.05  ~ "*",
                         TRUE ~ ""))

db_tests
len_tests

# y-positions for stars (above tallest bar at each bin)
len_ann <- len_sum %>% group_by(len) %>% summarise(y = max((mean + se)*100) + 1, .groups = "drop") %>%
  left_join(len_tests %>% select(len, sig), by = "len")
db_ann  <- db_sum  %>% group_by(db ) %>% summarise(y = max((mean + se)*100) + 1, .groups = "drop") %>%
  left_join(db_tests %>% select(db , sig), by = "db")

# --------------------------------------------
# Plots (style matches your sn-2 figure)
# --------------------------------------------
custom_greys <- c(groups_keep[1] = "#222222", groups_keep[2] = "#A9A9A9")
pd <- position_dodge(width = 0.9)

p_len_total <- ggplot(len_sum, aes(x = len, y = mean*100, fill = group)) +
  geom_col(position = pd, width = 0.85) +
  geom_errorbar(aes(ymin = (mean - se)*100, ymax = (mean + se)*100),
                position = pd, width = 0.85, linewidth = 0.3) +
  geom_vline(data = len_mark, aes(xintercept = mu, color = group),
             linetype = "dashed", linewidth = 0.6, alpha = 0.9, show.legend = FALSE) +
  #geom_text(data = len_ann, aes(x = len, y = y, label = sig), size = 3, vjust = 0) +
  scale_fill_manual(values = custom_greys, breaks = groups_keep, name = NULL) +
  scale_color_manual(values = custom_greys, guide = "none") +
  scale_x_continuous(limits = c(27,45), breaks = len_bins) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Total length (carbons)", y = "mol % phospholipids") +
  theme_pubr(base_size = 6) +
  theme(legend.position = "none")

p_db_total <- ggplot(db_sum, aes(x = db, y = mean*100, fill = group)) +
  geom_col(position = pd, width = 0.85) +
  geom_errorbar(aes(ymin = (mean - se)*100, ymax = (mean + se)*100),
                position = pd, width = 0.85, linewidth = 0.3) +
  geom_vline(data = db_mark, aes(xintercept = mu, color = group),
             linetype = "dashed", linewidth = 0.6, alpha = 0.9, show.legend = FALSE) +
  #geom_text(data = db_ann, aes(x = db, y = y, label = sig), size = 3, vjust = 0) +
  scale_fill_manual(values = custom_greys, breaks = groups_keep, name = NULL) +
  scale_color_manual(values = custom_greys, guide = "none") +
  scale_x_continuous(limits = c(-0.5,8.5), breaks = db_bins) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Total double bonds", y = "mol % phospholipids") +
  theme_pubr(base_size = 6) +
  theme(legend.position = "none")

p_len_total
ggsave(here("04-pdf", "k562_totallength.pdf"), width = 3, height = 3, units = "cm")

p_db_total
ggsave(here("04-pdf", "k562_totaldb.pdf"), width = 3, height = 3, units = "cm")

# OK, now let's do a fluidity analysis!
plfidata = pldata %>% 
  filter(!(class %in% class_minor)) %>% 
  calc_plfi()

# check that my weighting is correct
# yep looks fine!
plfidata %>% 
  group_by(eid) %>% 
  summarise(tm = weighted.mean(tm, frac_molar, na.rm = TRUE))

# visualize it like DBI
plfidata %>% 
  # average feature-wise within treatments
  filter(!is.na(plfi)) %>% 
  arrange(pres_gro, time_gro) %>% 
  mutate(
    trt = str_glue("{pres_gro} bar, {time_gro} h"),
    trt = trt %>% factor(levels = unique(.))
  ) %>% 
  group_by(id, tm, pres_gro, time_gro, trt, class, scheme) %>% 
  summarize(frac_molar = mean(frac_molar)) %>% 
  # and renorm
  group_by(pres_gro, time_gro, trt) %>% 
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  # and plot
  gg_plvisc(
    binwidth = 5,
    x = tm_bin,
    y = frac_molar,
    fill = class,
    rows = vars(trt)
  )

# plot PLFI mean/SEM just like PLCI
plfidata %>% 
  filter(!is.na(plfi)) %>% 
  mutate(class = class %>% factor(levels = ord_tm)) %>% 
  arrange(eid) %>% 
  group_by(eid, pres_gro, time_gro, class, scheme) %>% 
  summarise(plfi = sum(plfi)) %>% 
  # average by treatment and also get SEMs
  group_by(eid, pres_gro, time_gro, scheme) %>%
  mutate(plfi_tot = ifelse(row_number() == 1, sum(plfi), NA)) %>% 
  group_by(pres_gro, time_gro, class, scheme) %>%
  summarize(
    plfi = mean(plfi),
    plfi_sem = sd(plfi_tot, na.rm = TRUE) / sqrt(n()),
    plfi_tot = mean(plfi_tot, na.rm = TRUE)
  ) %>% 
  arrange(pres_gro, time_gro) %>% 
  mutate(
    trt = str_glue("{pres_gro} bar, {time_gro} h"),
    trt = trt %>% factor(levels = unique(.))
  ) %>% 
  gg_plfi(
    darkmode = FALSE,
    #rows = scheme,
    cols = trt,
    y = plfi
  ) +
  # there are two errorbar colors
  # so you can make the downward-pointing side white-on-black if you like
  geom_errorbar(
    aes(
      x = 0.3, # hardcoded into gg_plcurv()
      ymin = plfi_tot - plfi_sem,
      ymax = plfi_tot
    ),
    color = "black",
    #color = "white",
    width = 0
  ) +
  geom_errorbar(
    aes(
      x = 0.3, # hardcoded into gg_plcurv()
      ymin = plfi_tot,
      ymax = plfi_tot + plfi_sem
    ),
    color = "black",
    width = 0
  ) +
  theme_pubr() +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    strip.text.x = element_text(angle = 90),
    legend.position = "right"
  ) +
  labs(
    title = "K562 PLFI by treatment",
    fill = "Lipid class",
    x = "Treatment"
  )

ggsave(here("04-pdf", "k562_ms2_plfitrt_20250630a.pdf"), width = 4, height = 4)
ggsave(here("05-png", "k562_ms2_plfitrt_20250630a.png"), width = 4, height = 4)
