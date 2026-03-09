## Ingest c0 values from literature, generate estimates
## to apply to lipidomic data

library(tidyverse)
library(here)
library(httr)
library(ggpubr)

source(here("00-scripts", "seelipids_helpers.R")) 

file_curv = here("01-rawdata", "curvature.csv")
# this is a Google Sheets key for the curvature info file, to facilitate a nice Excel-like interface
gsht_curv = "1_LzSU0LEDQG5G3tDySbYWHwaKmEZoR2fGkye8CSc0QY"

# standard reference temperature for curvature contributions
temp_std = 20
temp_tol = 2.5

# refresh c0 spreadsheet from online interface
GET(str_glue("https://docs.google.com/spreadsheet/ccc?key={gsht_curv}&output=csv")) %>% 
  content("raw") %>% 
  write_file(file_curv)

# and load the data in
data_curv = file_curv %>% read_csv()

# let's just see all of them!
plot_allc0 = data_curv %>% 
  mutate(class = class %>% factor(levels = names(chroma_cl))) %>% 
  pivot_longer(
    cols = c("carbsn1", "dbonsn1"),
    names_to = "var",
    values_to = "val"
  ) %>% 
  arrange(var, class) %>% 
  mutate(varclass = paste(var, class) %>% factor(levels=unique(.))) %>%
  ggplot(
    aes(
      x = val,
      y = c0
    )
  ) +
  facet_wrap(~varclass, scale = "free_x", nrow = 2) +
  geom_point(
    aes(
      color = temp,
      shape = !is.na(relax)
    )
  ) +
  geom_smooth(method="lm", se=FALSE) +
  scale_colour_viridis_c() +
  scale_shape_manual(values = c(1, 16))

#plot_allc0

# filter to 0 bar values
curv_rated = data_curv %>% 
  filter(press == 0) %>% 
  # generate a quality score of sorts; higher is better
  mutate(
    temp_std = temp_std,
    goodness = 
      as.integer(!est) + # a point for not being estimated
      as.integer(!is.na(relax)) + # a point for having filler
      as.integer(str_detect(buff, "H2O")) + # a point for being in water
      as.integer(npln) + # one for being at neutral plane
      as.integer((temp == temp_std) | !is.na(tslope)) # a point for "standard temp" or at least a correction factor
  )

## SCHEME 1: 

# prioritize sn-1 18:1, try to adjust to 20°C
c0_oleoyl_best = curv_rated %>% 
  # simple acyl chain criteria
  filter(
    # sn-1 needs to be oleoyl
    (carbsn1 == 18) &
      (dbonsn1 == 1 ) |#&
      ## needs to be lyso-, or
      #(str_detect(class, 'L') |
      ## sn-2 also needs to be oleoyl
      #(carbsn2 == 18) &
      #(dbonsn2 == 1 )) |
      # dioleoyl plasmalogen not available!
      str_detect(class, "P-PE")
  ) %>% 
  arrange(class, -goodness) %>% 
  # grab just the top-quality data for each class
  group_by(class) %>% 
  filter(goodness == max(goodness)) %>% 
  ungroup() %>% 
  # and adjust for temperature
  mutate(
    temp_std = temp_std,
    c0  = c0 + ifelse(!is.na(tslope), tslope * (temp_std - temp), 0),
    tol = tol + abs(ifelse(!is.na(slopetol), ifelse(!is.na(tslope), tslope, 0) * (temp_std - temp), 0)),
    temp = ifelse(is.na(tslope), temp, temp_std)
  )
# then average
c0_oleoyl_mean = c0_oleoyl_best %>% 
  # average by class
  group_by(class) %>% 
  summarize(
    c0 = mean(c0),
    tol = max(tol, na.rm = TRUE),
    tol = ifelse(tol == -Inf, NA, tol),
    temp = mean(temp),
    across(contains("sn1"), mean),
    goodness = mean(goodness),
    ref = list(ref)
  )
# might want to add an adjustment here so P-PE looks "saturated"?
# or might not, for conservatism

## SCHEME 2

# sn-1 saturated, sn-2 oleoyl, try to adjust to 20°C.
c0_satsn1_best = curv_rated %>% 
  # simple acyl chain criteria
  filter(
    (dbonsn1 == 0) &
      ((carbsn2 == 18) & (dbonsn2 == 1) |
       str_detect(class, 'L')) |
      (class == "PS") # (no saturated PSs)
  ) %>% 
  arrange(class, -goodness) %>% 
  # grab just the top-quality data for each class
  group_by(class) %>% 
  filter(goodness == max(goodness)) %>% 
  ungroup() %>% 
  # and adjust for temperature
  mutate(
    temp_std = temp_std,
    c0  = c0 + ifelse(!is.na(tslope), tslope * (temp_std - temp), 0),
    tol = tol + abs(ifelse(!is.na(slopetol), ifelse(!is.na(tslope), tslope, 0) * (temp_std - temp), 0)),
    temp = ifelse(is.na(tslope), temp, temp_std)
  )
# then average
c0_satsn1_mean = c0_satsn1_best %>% 
  # average by class
  group_by(class) %>% 
  summarize(
    c0 = mean(c0),
    tol = max(tol, na.rm = TRUE),
    tol = ifelse(tol == -Inf, NA, tol),
    temp = mean(temp),
    across(contains("sn1"), mean),
    goodness = mean(goodness),
    ref = list(ref)
  )

## SCHEME 3

# build carbsn1 x dbonsn1 regressions for PE, PC only
# It's evident in the plot below that chain length and unsat effects
# on the lysolipids are not statistically significant,
# and just the means should be used.
# use only really good data
# classes where sn-1 chain will be considered
class_linreg = c("PC", "O-PC", "P-PC", "PE", "O-PE", "P-PE", "OO-PE", "PI") 
linreg_data = curv_rated %>% 
  filter(
    (class %in% class_linreg) &
      !is.na(relax) &
      str_detect(buff, "H2O") &
      (npln | str_detect(class, 'L')) & # none of the lysos are neutral-plane
      # also they were all measured at 22
      (between(temp, temp_std-temp_tol, temp_std+temp_tol) | !is.na(tslope))
      #((temp %in% c(20, 22)) | !is.na(tslope))
  ) %>% 
  # adjust for temperature
  mutate(
    temp_std = temp_std,
    # zero out the contribution if temp is >5°C off-target and cannot be adjusted
    #tslope = ifelse(is.na(tslope) & (abs(temp - temp_std) <= 5), 0, tslope),
    c0  = c0 + ifelse(!is.na(tslope), tslope * (temp_std - temp), 0),
    temp = ifelse(is.na(tslope), temp, temp_std)
  )

# let's see those groomed data!
plot_linreg = linreg_data %>% 
  ## just look at DO species
  #filter((carbsn1 == 18) & (dbonsn1 %in% c(0, 1))) %>% 
  pivot_longer(
    cols = c("carbsn1", "dbonsn1"),
    names_to = "var",
    values_to = "val"
  ) %>% 
  mutate(varclass = paste(var, class) %>% factor(levels=unique(.))) %>%
  ggplot(
    aes(
      x = val,
      y = c0,
      color = class,
      group = class
    )
  ) +
  facet_wrap(~var, scale = "free_x", nrow = 1) +
  geom_point() +
  geom_text(aes(label = ref)) +
  geom_smooth(method="lm", se=TRUE) +
  scale_color_manual(values = chroma_cl) +
  labs(
    title = "c0 by sn-1 chain length and unsat\nadjusted as necessary to 20°C"
  )
# by the looks of it, effect of +C2 is about the same as that of 1 dbond
# that means the coefficient should be about half.
plot_linreg

# a very simple model
# it would be nice if errors could be incorporated
linreg = linreg_data %>% 
  filter(class %in% c("PE", "PC")) %>% # leaving lysolipids out for reason given above
  #mutate(dbonsn1 = as.logical(dbonsn1)) %>% # code as categorical
  # improved, very simple!
  lm(c0 ~ class + dbonsn1 + carbsn1, .) 
  
summary(linreg)

# new class offsets
# for classes on which I want to impose sn-1 chain effects fitted from PE+PC
class_base = c("PE", "PC") # offsets are based on these
class_offs = c("O-PC", "P-PC", "O-PE", "P-PE", "OO-PE", "PI") # these classes get offsets
offsets = linreg_data %>% 
  # there is just one representative of each new class
  filter(class %in% class_offs) %>% 
  # find the most closely corresponding actual data
  left_join(
    linreg_data %>% 
      filter(class %in% class_base),
    by = c("carbsn1", "dbonsn1", "carbsn2", "dbonsn2"),
    relationship = "many-to-many",
    suffix = c("", "_base")
  ) %>% 
  # average the reference hits
  group_by(across(!contains("_base")), class_base) %>% 
  summarise(across(contains("_base") & where(is.numeric), mean)) %>% 
  # calculate the offset for each class
  group_by(class) %>% 
  mutate(offs = c0 - c0_base) %>% 
  # and select the offset from the class that is closest
  filter(abs(offs) == min(abs(offs))) %>% 
  # be sure to retain the reference class ID
  select(class, class_base, offs)

# Run out predictions of the base model (PE, PC only)
c0_chain_fx_base = crossing(
  class   = c("PC", "PE"),
  carbsn1 = seq(10, 22, 1),
  dbonsn1 = c(0, 1),
  temp    = 20
) %>% 
  rowwise() %>% 
  mutate(
    # run predictions w/90% CI
    c0  = predict(linreg, cur_data(), interval = "confidence", level = 0.90) %>% 
      as.list() %>% setNames(c("c0", "lo", "hi")) %>% list(),
    # insert the lm formula as ref
    ref = linreg$call %>% as.character() %>% .[[2]] %>% list()
  ) %>% 
  unnest_wider(c0) %>% 
  mutate(tol = c0-lo) %>% 
  select(-lo, -hi)

# Expand to extra classes using constant c0 offsets
c0_chain_fx_offs = c0_chain_fx_base %>% 
  rename(class_base = class) %>% 
  right_join(
    offsets, 
    by = "class_base",
    relationship = "many-to-many"
  ) %>% 
  arrange(class, carbsn1, dbonsn1) %>% 
  rowwise() %>% 
  mutate(
    c0 = c0 + offs,
    ref = c(unlist(ref), str_glue("offset from {class_base}")) %>% list()
  )

# combine the PE/PC and offset predictions
c0_chain_fx = bind_rows(
  c0_chain_fx_base,
  c0_chain_fx_offs
)

# plot those results
plot_curvmodel = c0_chain_fx %>% 
  ggplot(
    aes(
      x = carbsn1,
      y = c0,
      ymin = c0-tol,
      ymax = c0+tol,
      fill = class,
      color = class
    )
  ) +
  # perhaps should be point?
  #geom_col(position = position_dodge()) +
  geom_ribbon(color = NA, alpha = 0.2) +
  geom_line() +
  #geom_errorbar(position = position_dodge()) +
  facet_wrap(~dbonsn1) +
  #scale_fill_manual(values = chroma_cl) +
  scale_color_manual(values = chroma_cl) +
  scale_fill_manual(values = chroma_cl) +
  theme_pubr() +
  labs(
    title = "Predicted curvature as function of headgroup, sn-1 length and unsat."
  )
plot_curvmodel

# bind all the schemes together into one table
c0_allschemes = bind_rows(
  c0_oleoyl_mean %>% mutate(scheme = "oleoyl"),
  c0_satsn1_mean %>% mutate(scheme = "satsn1"),
  c0_chain_fx    %>% mutate(scheme = "linreg"),
) %>% 
  ## these are all curvatures at 0 bar! redundant here?
  #filter(press == 0) %>% 
  # fill in values for the linreg scheme with an average from other two schemes
  # this automatically picks them, since there is only one value per class
  # in the other schemes (PE, PC are the only ones w/variable sn-1 chains)
  complete(class, scheme) %>% 
  group_by(class) %>% 
  mutate(
    ref = ifelse(!is.na(c0), ref, cur_data() %>% .$ref %>% unlist() %>% unique() %>% unname()), 
    tol = ifelse(!is.na(c0), tol, cur_data() %>% .$tol %>% mean(na.rm=TRUE)), # will return NaN if all NA
    c0  = ifelse(!is.na(c0), c0,  cur_data() %>% .$c0  %>% mean(na.rm=TRUE)),
    #NTS 20230710 get model-based errors in here too!
  ) %>% 
  fill(carbsn1, dbonsn1, temp, .direction = "downup")

# plot representative curvatures for each class
# sn-1 chain struc fixed
c0_allschemes %>% 
  filter(
    (scheme == "linreg") &
      (carbsn1 == 18) &
      (dbonsn1 == 1)
  ) %>% 
  ungroup() %>% 
  arrange(c0) %>% 
  mutate(class = class %>% factor(., levels = .)) %>% 
  ggplot(
    aes(
      x = class,
      y = c0,
      fill = class
    )
  ) +
  geom_hline(color = "black", yintercept = 0) +
  geom_col() +
  theme_pubr() +
  theme(
    axis.line.y = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(values = chroma_cl) +
  scale_y_reverse() +
  coord_flip() +
  labs(
    y = "c0 (1/Å) for sn-1 C18:1",
    x = "Phospholipid class"
  )
ggsave(here("04-pdf", "c0_contribs_oleoyl_20250629a.pdf"), width = 4, height = 4)
ggsave(here("05-png", "c0_contribs_oleoyl_20250629a.png"), width = 4, height = 4)  
