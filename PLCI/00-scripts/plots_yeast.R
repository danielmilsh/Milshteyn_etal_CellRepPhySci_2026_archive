library(tidyverse)
library(here)
library(gridExtra)
library(httr)

source(here("00-scripts", "seelipids_helpers.R"))


data_raw1 = read_csv((here("01-rawdata", "yeast", "A1_Data_QS-22-816_orig.csv")))
data_raw2 = read_csv(here("01-rawdata", "yeast", "A3_Lipotype_Report_Budin (QS-21-585)_pmol.csv"))
data_raw <- full_join(data_raw1, data_raw2, by = c("feature", "class")) #combine data from 2 files
#data_raw = read_csv(file_data3)

data_transposed = data_raw %>% 
  mutate_all(as.character) %>% #make all characters
  rownames_to_column() %>% #switch rows and columns
  pivot_longer(cols = -rowname, names_to = 'variable', values_to = 'value') %>% 
  pivot_wider(values_from = value, names_from = rowname)
names(data_transposed) = head(data_transposed, 1) %>% as.character() #make first row as headers 
data_transposed = data_transposed %>% tail(nrow(.)-2) #remove row 1


data_long = data_transposed %>% 
  pivot_longer(cols = -feature, names_to = "cmpd", values_to = "pct") %>% 
  rename("sample" = feature) #%>%
#replace_na(list("pct" = 0)) %>% #replace all NA's w/ 0 
data_long$pct = as.character(data_long$pct)


#Cer 32:0;2
#separate by " ", left of " " is class, right is saturation
#in saturation, separate by "-", left is chain 1, right is chain 2, one digit right of ";" is isoform
#chain 1 length, chain 2 length, sum is total length
#chain 1 db, chain 2 db, sum is total db

data_long = data_long %>%
  #replace_na(list("pct" = 0)) %>%
  mutate(id = cmpd) %>%
  mutate(eid = sample)#make columns with full compound info by copy

data_long = data_long %>% separate(sample, into = c("strain", "sp", "depth", "replicate"), sep ='-') %>% #separate sample name info by "-
  separate(cmpd, into = c("class", "saturation"), sep =' ') %>%#separate compound infor into class and saturation
  separate(saturation, into = c("chain1", "chain2"), sep ='_') %>% #separate compound info into chain 1 and 2
  separate(chain1, into = c("chain1", "hydrsn1"), sep =';') %>% #get db info alone
  separate(chain2, into = c("chain2", "hydrsn2"), sep = ';') #%>% #same

#make total OH column
data_long$hydrsn1 = as.numeric(data_long$hydrsn1)
data_long$hydrsn2 = as.numeric(data_long$hydrsn2)
data_long = data_long %>% 
  mutate(total_OH = hydrsn1 + hydrsn2) # add column with sum of chain 1 and 2 OH

#make total DB column
data_long = data_long %>%
  mutate(chain1copy = chain1) %>%
  mutate(chain2copy = chain2) %>%
  separate(chain1copy, into = c("carbsn1", "dbonsn1"), sep = ':') %>%
  separate(chain2copy, into = c("carbsn2", "dbonsn2"), sep = ':') %>%
  #mutate(dbonsn2 = replace_na(dbonsn2, 0)) %>%
  mutate(dbonds = as.numeric(dbonsn1) + as.numeric(dbonsn2))

#make total length column
data_long = data_long %>%
  mutate(carbsn2 = replace_na(as.numeric(carbsn2), 0)) %>%
  mutate(carbon = as.numeric(carbsn1) + as.numeric(carbsn2)) %>%
  unite(chain, c(chain1, chain2), sep = "-", na.rm = TRUE, remove = FALSE)

#rearrange columns
data_long = data_long %>%
  select(id, everything()) %>%
  select(-pct, pct)

data_long <- data_long %>%
  mutate(depth = case_when(
    depth %in% c("1bar", "0bar") ~ 1,
    depth == "250bar" ~ 250,
    TRUE ~ as.numeric(NA)  # Handle unexpected values
  )) %>%
  mutate(strain = case_when(
    strain == "WT" ~ "W303",
    strain == "fad2" ~ "fad2", 
    strain == "FM628" ~ "FM628"
  ))

pl_major <- data_long %>%
  filter(class %in% c( "PE", "PC", "PA", "PS", "PI", "LPC", "LPE")) %>%
  group_by(strain, depth, replicate) %>%
  mutate(pct = as.numeric(pct)) %>%
  replace_na(list("pct" = 0)) %>%
  mutate(frac_molar = pct/sum(pct)) %>%
  # After calculating frac_molar per replicate, average it across replicates for each strain, depth, and id
  #group_by(strain, replicate, depth, id, class, carbon, dbonds, carbsn1, carbsn2, dbonsn1, dbonsn2) %>% # Group by class as well if you want to keep classes separate
  #summarise(frac_molar = mean(frac_molar)) %>% # Use summarise to collapse rows and remove 'replicate'
  #ungroup() %>% # Ungroup before the final normalization
  #filter(frac_molar >= 0.01) %>%# only those in abundance over 1%
  mutate(carbsn1 = as.numeric(carbsn1)) %>%              # Convert to numeric
  filter(!is.na(carbsn1) & carbsn1 %% 2 == 0) %>%        # Keep only even values
  group_by(strain, depth, replicate) %>% # Regroup to ensure final normalization sums to 1 per id, depth, strain
  mutate(frac_molar = frac_molar / sum(frac_molar)) %>%
  ungroup() %>%
  filter(strain %in% c("W303"))


pl_major <- pl_major %>%
  mutate(
    strain = factor(strain, levels = c("W303", "FM628")),
    class = factor(class, levels = c("LPC", "LPE" , "PS", "PC", "PA", "PI", "PE"))
  ) 

ggplot() +
  facet_wrap(~ replicate) +
  labs(x = "Pressure (bar)", y = "mole fraction") +
  scale_x_continuous(breaks = pl_major %>% pull(depth) %>% unique) +
  geom_col(
    data = pl_major, 
    aes(
      x = depth,
      y = frac_molar,
      fill = factor(class)
    )
  ) +
  scale_fill_manual(values = chroma_cl) #+
theme_pubclean()

pldata <- pl_major

pl_maj_summary <- pl_major %>% 
  #filter(strain %in% c("W303")) %>%
  filter(strain %in% c("W303", "FM628")) %>%
  #group_by(strain, depth, class) %>%
  #summarize(normed_frac = sum(frac_molar)) %>%
  group_by(strain, depth, class) %>% 
  summarize(
    avg_frac = mean(frac_molar),
    sdev_frac = sd(frac_molar),
    serr_frac = sdev_frac/sqrt(n())
  )

source(here("00-scripts", "c0op_helpers.R"))

# Thin out the PL data by removing any lipid class that is <0.5 mol% in all samples
class_minor = pldata %>% 
  #filter(str_detect(sp, "Ecol_wild")) %>% 
  group_by(class, eid, depth) %>% 
  summarise(frac_molar = sum(frac_molar)) %>% 
  filter(max(frac_molar) < 0.005) %>% 
  .$class %>% unique() %>% print()

# Estimate curvature contributions
# runs all 3 schemes by default
plcidata = pldata %>% 
  mutate(dbonsn1 = as.double(dbonsn1)) %>%
  mutate(carbsn1 = as.double(carbsn1)) %>%
  # remove chain data for classes that are not part of the linreg
  # mutate(
  #   carbsn1 = ifelse(!(class %in% class_linreg), NA, carbsn1),
  #   dbonsn1 = ifelse(!(class %in% class_linreg), NA, dbonsn1)
  # ) %>% 
  # duplicate the input data for each scheme to be evaluated
  # cross_join(
  #   tibble(scheme = scheme)
  # ) %>% 
  # join the headgroup-only schemes
  left_join(
    c0_allschemes %>% 
      filter(scheme %in% c("linreg")) %>% 
      select(scheme, class, carbsn1, dbonsn1, c0, tol),
    by = c("class", "carbsn1", "dbonsn1")
  ) %>%
  # mutate(
  #   c0  = ifelse(is.na(c0 ), c0.linreg, c0 ),
  #   tol = ifelse(is.na(tol), tol.linreg,  tol)
  # ) %>% 
  select(-contains(".linreg")) %>% 
  # finally, calculate the actual curvature contributions and error
  # `plci` stands for PhosphoLipid Curvature Index
  mutate(
    plci = c0  * frac_molar,
    ctol = tol * frac_molar # this is valid since there's no error on frac_molar
  ) %>% 
  # important!
  replace_na(list(plci = 0, ctol = 0)) %>% 
  #filter(plci != 0) %>% # dangerous, removes missed things from the legend
  mutate(class = factor(class, levels = class_order_c0))
  #mutate(carbsn1 = as.numeric(carbsn1)) %>%
  #mutate(dbonsn1 = as.numeric(dbonsn1)) %>%
  #filter(str_detect(sp, "Ecol_wild")) %>% 
  #filter(!(class %in% class_minor)) %>% 
  #calc_plci()



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
    title = "Yeast PL curvature",
    fill = "Lipid class",
    x = "Treatment"
  )


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
    strip.text.x = element_text(angle = 90),
    legend.position = "right"
  ) +
  labs(
    title = "Ecoli PL curvature by treatment",
    fill = "Lipid class",
    x = "Incubation pressure"
  )
