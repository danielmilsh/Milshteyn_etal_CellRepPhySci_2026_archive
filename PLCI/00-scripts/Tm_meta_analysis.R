## Ingest Tm values from literature, generate estimates
## to apply to lipidomic data

library(tidyverse)
library(here)
library(httr)
library(broom)
library(ggpubr)

source(here("00-scripts", "seelipids_helpers.R")) 

# helper
tidyglance = function(mod){
  bind_cols(tidy(mod), glance(mod) %>% select(-statistic, -p.value))
}

file_melt = here("01-rawdata", "melttemps.csv")
# this is a Google Sheets key for the meltature info file, to facilitate a nice Excel-like interface
gsht_melt = "1sNhLxS7h5s99fw18LVCfE3cg-U6VxboRj4qGDIJyoTw"

# refresh c0 spreadsheet from online interface
GET(str_glue("https://docs.google.com/spreadsheet/ccc?key={gsht_melt}&output=csv")) %>% 
  content("raw") %>% 
  write_file(file_melt)

# and load the data in
data_melt = file_melt %>% read_csv()

# Simple fully additive model, 5 parameters
# the question is whether to treat the double bonds per chain
# as continuous or categorical, since they have nonmonotonic effects
modcompare_plfi = data_melt %>%
  # offset makes DSPX lipid the reference point
  #mutate(across(contains("carb"), ~.x-18)) %>% 
  #mutate(class = class %>% factor(levels = c("PC", "PG", "PS", "CL", "PE", "PA"))) %>% 
  {bind_rows(
    tibble(
      scheme = "monotonic",
      mod = {.} %>% 
        lm(Tm ~ class + carbsn1 + dbonsn1 + carbsn2 + dbonsn2 + 0, .) %>% list()
    ),
    tibble(
      scheme = "sn2_nmono",
      mod = {.} %>% 
        mutate(dbonsn2 = factor(dbonsn2)) %>%
        lm(Tm ~ class + carbsn1 + dbonsn1 + carbsn2 + dbonsn2 + 0, .) %>% list()
    ),
    tibble(
      scheme = "bothnmono",
      mod = {.} %>% 
        mutate(across(contains("dbon"), factor)) %>%
        lm(Tm ~ class + carbsn1 + dbonsn1 + carbsn2 + dbonsn2 + 0, .) %>% list()
    ),
    tibble(
      scheme = "logmonotonic",
      mod = {.} %>% 
        lm(Tm ~ class + log(carbsn1) + dbonsn1 + log(carbsn2) + dbonsn2 + 0, .) %>% list()
    ),
    tibble(
      scheme = "logsn2_nmono",
      mod = {.} %>% 
        mutate(dbonsn2 = factor(dbonsn2)) %>%
        lm(Tm ~ class + log(carbsn1) + dbonsn1 + log(carbsn2) + dbonsn2 + 0, .) %>% list()
    ),
    tibble(
      scheme = "logbothnmono",
      mod = {.} %>% 
        mutate(across(contains("dbon"), factor)) %>%
        lm(Tm ~ class + log(carbsn1) + dbonsn1 + log(carbsn2) + dbonsn2 + 0, .) %>% list()
    ),
    tibble(
      scheme = "logcanonical",
      mod = {.} %>% 
        # only take omega-9 and omega-3 sn-2 UFAs
        filter(
          (dbonsn2 == 0) | 
            ((dbonsn2 == 1) & (omegsn2 == 9)) |
            ((dbonsn2 %in% c(3, 5, 6)) & (omegsn2 == 3)) |
            ((dbonsn2 %in% c(2, 4)) & (omegsn2 == 6))
        ) %>% 
        mutate(across(contains("dbon"), factor)) %>%
        lm(Tm ~ class + log(carbsn1) + dbonsn1 + log(carbsn2) + dbonsn2 + 0, .) %>% list()
    ),
    # what about intersections?
    tibble(
      scheme = "logintersect",
      mod = {.} %>% 
        mutate(across(contains("dbon"), factor)) %>%
        lm(Tm ~ class + log(carbsn1) * dbonsn1 + log(carbsn2) * dbonsn2 + 0, .) %>% list()
    ),
  )} %>% 
  rowwise() %>% 
  mutate(tg = tidyglance(mod) %>% list()) %>% 
  unnest(tg)

# pull out just the models and their GoF measures
mods_plfi = modcompare_plfi %>% 
  group_by(scheme) %>% 
  select(scheme, mod, adj.r.squared, logLik, AIC, BIC, nobs) %>% 
  slice(1) %>% 
  arrange(AIC)

# look at the coeffs
modcompare_plfi %>% 
  filter(scheme == "logcanonical") %>% View()

# do some visual sanity checks
scheme_check = "logintersect"
#scheme_check = mods_plfi %>% ungroup() %>% slice(1) %>% .$scheme

# this is the disaturated plot on the Avanti page
tibble(
  class = c("PC", "PE", "PS", "PG"),
) %>% 
  cross_join(
    tibble(
      carbsn1 = seq(12, 24, 2),
      dbonsn1 = 0,
      carbsn2 = seq(12, 24, 2),
      dbonsn2 = 0
    )
  ) %>% 
  # necessary data prep
  #mutate(across(contains("carb"), ~.x-18)) %>% 
  mutate(across(contains("dbon"), factor)) %>%
  # use the best model yet
  cross_join(mods_plfi %>% filter(scheme == scheme_check)) %>% 
  ungroup() %>% 
  mutate(Tm_pred = predict(mod[[1]], newdata = cur_data())) %>% 
  ggplot(
    aes(
      x = carbsn1,
      y = Tm_pred,
      color = class,
      group = class
    )
  ) +
  geom_line() +
  theme_pubr() +
  scale_color_manual(values = chroma_cl) +
  labs(
    title = str_glue("Model: {scheme_check}"),
    x = "sn-1, 2 chain length (saturated)",
    y = "Predicted Tm (°C)"
  )

# This is the stearoyl-PUFA plot in Niebylski and Salem (Fig. 2)
crossing(
  class = c("PC"),
  carbsn1 = 18,
  dbonsn1 = 0,
) %>% 
  cross_join(
    bind_rows(
      crossing(
        carbsn2 = 18,
        dbonsn2 = seq(1, 3)
      ),
      crossing(
        carbsn2 = 20,
        dbonsn2 = seq(2, 5)
      ),
      crossing(
        carbsn2 = 22,
        dbonsn2 = seq(4, 6)
      )
    )
  ) %>% 
  # necessary data prep
  #mutate(across(contains("carb"), ~.x-18)) %>% 
  mutate(across(contains("dbon"), factor)) %>%
  # use the best model yet
  cross_join(
    mods_plfi %>% 
      filter(!str_detect(scheme, "mono") | str_detect(scheme, "both"))
  ) %>% 
  group_by(scheme) %>% 
  mutate(Tm_pred = predict(mod[[1]], newdata = cur_data())) %>% 
  ggplot(
    aes(
      x = dbonsn2,
      y = Tm_pred,
      color = class,
      group = paste(class, carbsn2)
    )
  ) +
  facet_wrap(~scheme, nrow = 1) +
  geom_line() +
  theme_pubr() +
  scale_color_manual(values = chroma_cl) +
  labs(
    title = "Compare Niebylski and Salem",
    x = "sn-2 double bonds; sn-1 18:0",
    y = "Predicted Tm (°C)"
  )

# logcanonical and bothnmono both darn good!

#Visually comparing to Niebylski and Salem Fig. 2, `logintersect` appears to be the most accurate, _however_ it is severely rank-deficient for sn-1 PUFAs (surprise). The interaction effects that are present are non-significant.
#
#Tm of SOPC appears high compared to Niebylski and Salem bc the value from NIST/Avanti is much higher (6°C vs. 0°C). I trust high values more; depression could be from impurity.