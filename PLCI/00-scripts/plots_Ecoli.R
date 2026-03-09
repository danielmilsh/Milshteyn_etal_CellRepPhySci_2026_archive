library(tidyverse)
library(here)
library(gridExtra)
library(httr)

source(here("00-scripts", "seelipids_helpers.R"))
source(here("00-scripts", "parse_lipidmaps_Ecoli.R"))
source(here("00-scripts", "c0op_helpers.R"))

# load the parsed data
file_tidy = here("02-tidydata", "Ecoli_lipidmaps_tidy.tsv")
file_tidy = here("02-tidydata", "lipidomics_data.csv")
file_tidy = here("02-tidydata", "MG_lipidomics_data_feb.csv")

# skip if parsing was just done and lmapdata_long is already in the global environment
lmapdata_long = read_csv(file_tidy)

# this file annotates unique extract IDs (`eid`) with metadata (`sp`, `treatment`, etc)
file_meta = here("01-rawdata", "Ecoli" , "Ecoli_lipidmaps_metadata.tsv")

file_meta = here("01-rawdata", "Ecoli" , "MG1566meta.csv")

# load metadata
metadata = file_meta %>% read_csv()

# join metadata to lipid data
lmapdata = lmapdata_long %>%
  ungroup() %>% 
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
  arrange(eid, id)

# store it
lmapdata %>% write_tsv(here("02-tidydata", "MG_lipidmaps_wmeta.tsv"))

# # where to store metadata
# #file_metdat = here("01-metadata", "metadata.csv")
# 
# # load the parsed data
# #file_tidy = here("02-tidydata", "lipidmaps_tidy.tsv")
# # skip if parsing was just done and lmapdata_long is already in the global environment
# #lmapdata_abs = read_tsv(file_tidy)
# 
# # this is a Google Sheets key for the metadata file, to facilitate a nice Excel-like interface
# #gsht_metdat = "1w06ib3W2fOUGj7iOjV9cvEZ4mWMnG1YH1deafM71KgQ"
# 
# # refresh c0 spreadsheet from online interface
# #  GET(str_glue("https://docs.google.com/spreadsheet/ccc?key={gsht_metdat}&output=csv")) %>%
# #  content("raw") %>%
# #  write_file(file_metdat)
# 
# # and load the data in
# #metadata = file_metdat %>% read_csv()
# 
# # join metadata to lipid data
# lmapdata_raw = lmapdata_abs %>%
#   left_join(metadata, by = "eid") %>%
#   # then make compound ID a factor so it plots in a uniform order
#   mutate(
#     # put headgroups in the order specified in the color palette,
#     # which can be changed in seelipids_helpers.R
#     class = factor(class, levels = names(chroma_cl)),
#     # then order things by total chain length and unsaturation
#     id = factor(
#       id,
#       levels = cur_data() %>%
#         distinct(id, class, carbon, dbonds) %>%
#         arrange(class, carbon, dbonds) %>%
#         pull(id) %>%
#         unique()
#     )
#   ) %>%
#   # sort the rows for easy viewing
#   arrange(depth, rep) %>%
#   mutate(eid = eid %>% factor(levels = unique(eid))) %>%
#   arrange(eid, id)
# 
# ## we have 39:0 and 39:1 PGs at fairly high abundance
# ## but I am gonna leave them alone for now
# #unfragd = lmapdata_raw %>%
# #  filter(id == "PE 41:6_19.26") %>%
# #  select(gn, sp, eid, class, frac_molar, dbonds, carbon, temp_stor) %>%
# #  cross_join(
# #    tibble(
# #      id = "PE 19:0/22:6_19.26",
# #      annot = "19:0/22:6",
# #      carbsn1 = 19,
# #      carbsn2 = 22,
# #      dbonsn1 = 0,
# #      dbonsn2 = 6
# #    )
# #  )
# 
# # other than that we are only considering species with MS2-confirmed stereochem.
# lmapdata = lmapdata_raw %>%
#   filter(ms2)# %>%
#   #bind_rows(unfragd)
# 
# # store it
# lmapdata %>% write_tsv(here("02-tidydata", "lipidmaps_wmeta.tsv"))

#### PLOTS

library(dplyr)
library(stringr)
library(dplyr)

# how many eid groups exist total?
n_eid <- n_distinct(lmapdata$eid)

pldata <- lmapdata %>%
  filter(class %in% c("PE", "PG", "CL")) %>% #,
         #carbsn1 > 0) %>%
  
  # count how many eids each lipid id appears in
  group_by(id) %>%
  filter(n_distinct(eid) == n_eid) %>%   # <-- KEY STEP
  ungroup() %>%
  
  # renormalize AFTER filtering
  group_by(eid) %>%
  mutate(frac_molar = rab / sum(rab)) %>%
  ungroup()

write.csv(pldata, "C:/Users/danie/OneDrive/Desktop/Manuscripts/YeastHomeocurvature/PNAS/CellRepPhysSci_revision/MG1655lipidomics.csv")

  # confirm each eid has the same lipid set
pldata_common %>% count(eid)  # rows per eid (should be identical if other cols match)
pldata_common %>% group_by(eid) %>% summarise(sum_frac = sum(frac_molar))

pldata %>% ungroup() %>% distinct(class, annot) %>% nrow()

pldata %>% 
  #filter(str_detect(sp, "Ecol_wild")) %>% 
  gg_headgp(
    darkmode = FALSE,
    # aesthetic mappings are passed straight thru
    x = eid,
    y = frac_molar,
    fill = class
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(
    title = "Ecoli phospholipids",
    x = "Sample ID",
    y = "Mole fraction phospholipids",
    fill = "Lipid class"
  )
# save vector and raster images
ggsave(here("04-pdf", "MG1566_phospholipids_ms2.pdf"), width = 8, height = 5)
ggsave(here("05-png", "Ecoli_phospholipids_ms2_20250618a.png"), width = 8, height = 5)

class_foci = c("PE", "PG", "CL")

# acyl chain comparisons
pldata %>% 
  #filter(str_detect(sp, "Ecol_wild")) %>% 
  filter(class %in% class_foci) %>% 
  arrange(eid) %>% 
  # now add all the lipids in every class and #C together
  group_by(sp, eid, depth, class, carbon) %>% 
  summarize(frac_molar = sum(frac_molar)) %>% 
  # now average by genus
  group_by(sp, depth, class, carbon) %>% 
  summarize(frac_molar = mean(frac_molar)) %>%
  # normalize within each headgroup
  group_by(sp, depth, class) %>%
  mutate(frac_molar = frac_molar/max(frac_molar)) %>% 
  ungroup() %>% 
  arrange(depth) %>% 
  mutate(
    trt = str_glue("{depth} bar"),
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
    title = "Ecoli chain lengths",
    x = "Total acyl carbons",
    y = "Relative frequency"
  )
# save vector and raster images
ggsave(here("04-pdf", "Ecoli_ms2_acylcarbons_20250625a.pdf"), width = 8, height = 4)
ggsave(here("05-png", "Ecoli_ms2_acylcarbons_20250625a.png"), width = 8, height = 4)

pldata %>% 
  #filter(str_detect(sp, "Ecol_wild")) %>% 
  filter(class %in% class_foci) %>% 
  arrange(eid) %>% 
  # now add all the lipids in every class and #C together
  group_by(sp, eid, depth, class, dbonds) %>% 
  summarize(frac_molar = sum(frac_molar)) %>% 
  # now average by genus
  group_by(sp, depth, class, dbonds) %>% 
  summarize(frac_molar = mean(frac_molar)) %>%
  # normalize within each headgroup
  group_by(sp, depth, class) %>%
  mutate(frac_molar = frac_molar/max(frac_molar)) %>% 
  ungroup() %>% 
  arrange(depth) %>% 
  mutate(
    trt = str_glue("{depth} bar"),
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
    rows = vars(depth)
  ) +
  #coord_flip() +
  labs(
    title = "Ecoli unsaturations",
    x = "Total unsaturations",
    y = "Relative frequency"
  )
# save vector and raster images
ggsave(here("04-pdf", "Ecoli_ms2_acyldbonds_20250625a.pdf"), width = 8, height = 4)
ggsave(here("05-png", "Ecoli_ms2_acyldbonds_20250625a.png"), width = 8, height = 4)

# Thin out the PL data by removing any lipid class that is <0.5 mol% in all samples
class_minor = pldata %>% 
  #filter(str_detect(sp, "Ecol_wild")) %>% 
  group_by(class, eid, depth) %>% 
  summarise(frac_molar = sum(frac_molar)) %>% 
  #filter(max(frac_molar) < 0.005) %>% 
  .$class %>% unique() %>% print()

# Estimate curvature contributions
# runs all 3 schemes by default
plcidata = pldata %>% 
  #filter(str_detect(sp, "Ecol_wild")) %>% 
  #filter(!(class %in% class_minor)) %>% 
  calc_plci()

# plot PLCI by individual
plcidata %>% 
  filter(scheme == "linreg") %>% 
  arrange(eid) %>% 
  group_by(eid, depth, class, scheme) %>% #filter(class=="LPC") %>% View()
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
    title = "Ecoli PL curvature by sample",
    fill = "Lipid class",
    x = "Sample"
  )

ggsave(here("04-pdf", "MGEcoli_ms2_plcisamp_20250626a.pdf"), width = 6, height = 4)
ggsave(here("05-png", "Ecoli_ms2_plcisamp_20250626a.png"), width = 6, height = 4)

# plot PLCI mean/SEM
plcidata %>% 
  filter(scheme == "linreg") %>% 
  arrange(eid) %>% 
  group_by(eid, depth, class, scheme) %>% 
  summarise(plci = sum(plci)) %>% 
  # average by treatment and also get SEMs
  group_by(eid, depth, scheme) %>%
  mutate(plci_tot = ifelse(row_number() == 1, sum(plci), NA)) %>% 
  group_by(depth, class, scheme) %>%
  summarize(
    plci = mean(plci),
    plci_sem = sd(plci_tot, na.rm = TRUE) / sqrt(n()),
    plci_tot = mean(plci_tot, na.rm = TRUE)
  ) %>% 
  arrange(depth) %>% 
  mutate(
    trt = str_glue("{depth} bar"),
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
    title = "Ecoli PL curvature",
    fill = "Lipid class",
    x = "Treatment"
  )

ggsave(here("04-pdf", "Ecoli_plcitrt_20250626b.pdf"), width = 4, height = 4)
ggsave(here("05-png", "Ecoli_plcitrt_20250626b.png"), width = 4, height = 4)

# t-tests
plci_tot = plcidata %>% 
  filter(scheme == "linreg") %>% 
  arrange(eid) %>% 
  # calc total PLCI
  group_by(eid, depth, scheme) %>% 
  summarise(plci = sum(plci)) %>% 
  ungroup()

# 6 h @ 1 vs. 250 bar
plci_tot %>% 
  #filter(time_gro == 6) %>% 
  pivot_wider(names_from = "depth", values_from = "plci") %>% 
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
  #filter(time_gro == 6) %>% 
  pivot_wider(names_from = "depth", values_from = "plci") %>% 
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

# OK, try plotting the pooled data
plcidata %>% 
  filter(scheme == "linreg") %>% 
  arrange(eid) %>% 
  group_by(eid, depth, class, scheme) %>% 
  summarise(plci = sum(plci)) %>% 
  # average by treatment and also get SEMs
  group_by(eid, depth, scheme) %>%
  mutate(plci_tot = ifelse(row_number() == 1, sum(plci), NA)) %>% 
  group_by(depth, class, scheme) %>%
  summarize(
    plci = mean(plci),
    plci_sem = sd(plci_tot, na.rm = TRUE) / sqrt(n()),
    plci_tot = mean(plci_tot, na.rm = TRUE)
  ) %>% 
  arrange(depth) %>% 
  mutate(
    trt = str_glue("{depth} bar"),
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
    strip.text.x = element_text(angle = 0),
    legend.position = "right"
  ) +
  labs(
    #title = "Ecoli PL curvature by treatment",
    fill = "Lipid class",
    x = "Incubation pressure"
  )

ggsave(here("04-pdf", "MGEcoli_ms2_plcipress.pdf"), width = 4, height = 4)
ggsave(here("05-png", "Ecoli_ms2_plcipress_20250626a.png"), width = 4, height = 4)



library(paletteer)
# Plot with the customized facet order
classes <- plcidata %>% 
  #filter(time_gro == 6) %>%
  filter(scheme == "linreg") %>% 
  arrange(eid) %>% 
  group_by(eid, depth, class, scheme) %>% 
  summarise(avg_molar = sum(frac_molar)) %>% 
  # average by treatment and also get SEMs
  #group_by(eid, pres_gro, time_gro, scheme) %>%
  #mutate(frac_tot = ifelse(row_number() == 1, sum(frac_molar), NA)) %>% 
  group_by(depth, class, scheme) %>%
  summarize(
    frac_molar = mean(avg_molar),
    sem = sd(avg_molar) / sqrt(n()),
    tot = mean(avg_molar, na.rm = TRUE)
  ) #%>%
ggplot() +
  #facet_wrap(~ pres_gro) +
  labs(x = "Pressure (bar)", y = "mole fraction") +
  scale_x_continuous(breaks = pldata %>% pull(depth) %>% unique) +
  geom_col(
    data = classes, 
    aes(
      x = depth,
      y = frac_molar,
      #fill = class
      fill = factor(class, levels = c("LPC", "PS", "LPE", "PG", "PC", "O-PC", "P-PC", "PA", "PI", "PE", "P-PE", "OO-PE", "CL"))
    )
  ) +
  labs(fill = "Lipid") +
  theme_pubclean(base_size = 10) +
  scale_fill_manual(values = chroma_cl) +  
  #scale_fill_paletteer_d(("ggsci::category20_d3")) + #scale_fill_manual(values = chroma_cl) +
  #scale_fill_locuszoom() +
  theme(legend.position = "right") #+ 
theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())
ggsave(here("04-pdf", "MGEcoli_classcomp_20250726a.pdf"), width = 2, height = 4)



#meow
class_levels <- c("LPC","PS","LPE","PG","PC","O-PC","P-PC","PA","PI","PE","P-PE","OO-PE","CL")

classes2 <- classes %>%
  mutate(
    class = factor(class, levels = class_levels),
    depth = factor(depth, levels = sort(unique(depth)))
  )

ggplot(classes2, aes(x = depth, y = mean, fill = class)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = pmax(mean - sem, 0), ymax = mean + sem),
                width = 0.2, linewidth = 0.25) +
  facet_wrap(~ class, scales = "free_y") +
  scale_fill_manual(values = chroma_cl, guide = "none") +
  labs(x = "Depth", y = "mole fraction") +
  theme_pubclean(base_size = 6)

######################################################


#total db comp

db_maj_summary <- plcidata %>% 
  filter(scheme == "linreg") %>%
  #filter(carbsn1 > 0)  %>%
  filter(dbonds < 3) %>%
  arrange(eid) %>%
  group_by(eid, depth, dbonds) %>%
  summarize(db_frac = sum(frac_molar)) %>%
  group_by(depth, dbonds) %>% 
  summarize(
    avg_frac = mean(db_frac),
    sdev_frac = sd(db_frac),
    serr_frac = sdev_frac/sqrt(n()),
  )


dodge <- position_dodge(width = 0.9)  # should match your column dodging

db_maj_summary %>% 
  ggplot(aes(x = dbonds, y = avg_frac*100, fill = as.factor(depth))) +
  geom_col(position = dodge) +
  geom_errorbar(
    aes(
      ymin = (avg_frac - serr_frac)*100,
      ymax = (avg_frac + serr_frac)*100
    ),
    position = dodge, width = 0,
    width = 0.4
  ) +
  labs(
    x = "Total Double Bonds",
    y = "mol % phospolipids", 
    fill = "Pressure (bar)"
  ) +
  theme_pubr(base_size = 8) +
  scale_fill_grey() +
  theme(legend.position = "none")

ggsave(here("04-pdf", "MGEcoli_dbcomp.pdf"), width = 4, height = 4)



db_rep <- plcidata %>% 
  filter(scheme == "linreg",
         dbonds < 3) %>%
  group_by(eid, depth, dbonds) %>%
  summarize(db_frac = sum(frac_molar), .groups = "drop")
library(dplyr)
library(broom)

ttest_results <- db_rep %>%
  group_by(dbonds) %>%
  do(tidy(t.test(db_frac ~ depth, data = .))) %>%
  ungroup()
ttest_results
ttest_results <- ttest_results %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"))





####################################################


#total length comp
len_maj_summary <- plcidata %>% 
  filter(scheme == "linreg") %>%
  #filter(time_gro == 6)  %>%
  #filter(total_DB >= 0 & total_DB <= 4) %>%
  arrange(eid) %>%
  group_by(eid, depth, carbon, scheme) %>%
  summarize(len_frac = sum(frac_molar)) %>%
  group_by(depth, carbon) %>% 
  summarize(
    avg_frac = mean(len_frac),
    sdev_frac = sd(len_frac),
    serr_frac = sdev_frac/sqrt(n()),
  )

library(ggbreak) 
library(patchwork)

len_maj_summary %>% 
  #filter(time_gro == 6) %>%
  filter(carbon < 39) %>% ###double check filtering is with %in%
  #filter(carbon %% 2 == 0) %>%
  ggplot() +
  #facet_wrap(~eid) + 
  geom_col(
    #data = pl_maj_summary, 
    position = position_dodge(),
    aes(
      x = carbon,
      y = avg_frac*100,
      fill = as.factor(depth)
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
    position = position_dodge(0.9), width = 0,
    aes(
      x = carbon,
      ymin = (avg_frac - serr_frac)*100,
      ymax = (avg_frac + serr_frac)*100,
      group = depth
    )
  ) +
  theme(legend.position = 'none') +
  scale_x_continuous(limits = c(27, 39), breaks = seq(28, 38, 2)) #+  # Adjusted to fit your OD values
  #scale_x_break(c(0, ), ticklabels = c(28, 30, 32, 34, 36, 38)) +
  #scale_x_break(c(38, 64), ticklabels = c(28, 30, 32, 34, 36, 38, 62, 64, 66, 68, 70))
ggsave(here("04-pdf", "MGEcoli_ms2_lengthcomp.pdf"), width = 4, height = 4)


len_rep <- plcidata %>% 
  filter(scheme == "linreg",
         carbon < 39) %>%
  group_by(eid, depth, carbon) %>%
  summarize(len_frac = sum(frac_molar), .groups = "drop")
library(dplyr)
library(broom)

ttest_results <- len_rep %>%
  group_by(carbon) %>%
  do(tidy(t.test(len_frac ~ depth, data = .))) %>%
  ungroup()
ttest_results
ttest_results <- ttest_results %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"))



##########################################################################

# OK, now let's do a fluidity analysis!
plfidata = pldata %>% 
  #filter(!(class %in% class_minor)) %>% 
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
  arrange(depth) %>% 
  mutate(
    trt = str_glue("{depth} bar"),
    trt = trt %>% factor(levels = unique(.))
  ) %>% 
  group_by(id, tm, depth, trt, class, scheme) %>% 
  summarize(frac_molar = mean(frac_molar)) %>% 
  # and renorm
  group_by(depth, trt) %>% 
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
  group_by(eid, depth, class, scheme) %>% 
  summarise(plfi = sum(plfi)) %>% 
  # average by treatment and also get SEMs
  group_by(eid, depth, scheme) %>%
  mutate(plfi_tot = ifelse(row_number() == 1, sum(plfi), NA)) %>% 
  group_by(depth, class, scheme) %>%
  summarize(
    plfi = mean(plfi),
    plfi_sem = sd(plfi_tot, na.rm = TRUE) / sqrt(n()),
    plfi_tot = mean(plfi_tot, na.rm = TRUE)
  ) %>% 
  arrange(depth) %>% 
  mutate(
    trt = str_glue("{depth} bar"),
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
    title = "Ecoli PLFI by treatment",
    fill = "Lipid class",
    x = "Treatment"
  )

ggsave(here("04-pdf", "Ecoli_ms2_plfitrt_20250630a.pdf"), width = 4, height = 4)
ggsave(here("05-png", "Ecoli_ms2_plfitrt_20250630a.png"), width = 4, height = 4)
