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
  mutate(across(contains("carb"), ~.x-18)) %>% 
  {bind_rows(
    {.} %>% 
      lm(Tm ~ class + carbsn1 + dbonsn1 + carbsn2 + dbonsn2 + 0, .) %>% 
      tidyglance() %>% mutate(mod = "monotonic"),
    {.} %>% 
      mutate(dbonsn2 = factor(dbonsn2)) %>% 
      lm(Tm ~ class + carbsn1 + dbonsn1 + carbsn2 + dbonsn2 + 0, .) %>% 
      tidyglance() %>% mutate(mod = "sn2_nmono"),
    {.} %>% 
      mutate(across(contains("dbon"), factor)) %>% 
      lm(Tm ~ class + carbsn1 + dbonsn1 + carbsn2 + dbonsn2 + 0, .) %>% 
      tidyglance() %>% mutate(mod = "bothnmono"),
  )}

modcompare_plfi %>% 
  select(mod, adj.r.squared, logLik, AIC, BIC) %>% 
  distinct()
