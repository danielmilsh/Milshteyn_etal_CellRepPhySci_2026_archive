library(readr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggsci)
library(gridExtra)
library(tidyr)
library(multcomp)
library(here)


file_data1 = "~/yeast_lipotype032023/A1_Data_QS-22-816_orig.csv" #W303
file_data2 = "~/yeast_lipotype032023/A3_Lipotype_Report_Budin (QS-21-585)_pmol.csv" #FM628
file_data3 = "~/A1_Data_Q-24-1083.csv" #PEPC_yeast
file_data4 = "~/W303_FM628_PL_molpct.csv" #molpct W303 and FM628 combined
data_raw1 = read_csv(file_data1)
data_raw2 = read_csv(file_data2)
data_raw <- full_join(data_raw1, data_raw2, by = c("feature", "class")) #combine data from 2 files
data_raw = read_csv(file_data3)

source("~/seelipids_helpers.R")

##########################################################################################################################
####################################################  curate data  #######################################################
##########################################################################################################################

#run lines 1-115 for data parsing and PL summary

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
  mutate(lipid = cmpd) %>%
  mutate(full_id = sample)#make columns with full compound info by copy
  
data_long = data_long %>% separate(sample, into = c("strain", "species", "pressure", "replicate"), sep ='-') %>% #separate sample name info by "-
  separate(cmpd, into = c("class", "saturation"), sep =' ') %>%#separate compound infor into class and saturation
  separate(saturation, into = c("chain1", "chain2"), sep ='_') %>% #separate compound info into chain 1 and 2
  separate(chain1, into = c("chain1", "chain1_OH"), sep =';') %>% #get db info alone
  separate(chain2, into = c("chain2", "chain2_OH"), sep = ';') #%>% #same
#make pressure column only numerical
data_long$pressure = parse_number(as.character(data_long$pressure))
data_long = data_long %>%   # Replacing values
  mutate(pressure = replace(pressure, pressure == 0, 1))
#make total OH column
data_long$chain1_OH = as.numeric(data_long$chain1_OH)
data_long$chain2_OH = as.numeric(data_long$chain2_OH)
data_long = data_long %>% 
  mutate(total_OH = chain1_OH + chain2_OH) # add column with sum of chain 1 and 2 OH

#make total DB column
data_long = data_long %>%
  mutate(chain1copy = chain1) %>%
  mutate(chain2copy = chain2) %>%
  separate(chain1copy, into = c("chain1_length", "chain1_db"), sep = ':') %>%
  separate(chain2copy, into = c("chain2_length", "chain2_db"), sep = ':') %>%
  #mutate(chain2_db = replace_na(chain2_db, 0)) %>%
  mutate(total_DB = as.numeric(chain1_db) + as.numeric(chain2_db))

#make total length column
data_long = data_long %>%
  mutate(chain2_length = replace_na(as.numeric(chain2_length), 0)) %>%
  mutate(total_length = as.numeric(chain1_length) + as.numeric(chain2_length)) %>%
  unite(chain, c(chain1, chain2), sep = "-", na.rm = TRUE, remove = FALSE)

#rearrange columns
data_long = data_long %>%
  select(lipid, everything()) %>%
  select(-pct, pct)



##########################################################################################################################
##################################################  normalize data #######################################################
##########################################################################################################################
#major lipids filter
#filter phospholipids above threshold 
pl_major <- data_long %>%
  #filter(avg_pct >= 1) %>% # only those in abundance over 1%
  filter(class %in% c( "PE", "PC", "PA", "PS", "PI", "LPE", "LPC", "LPI", "LPS")) %>%
  #filter(strain %in% c("WT", "FM628")) %>%
  group_by(strain, pressure, replicate) %>% 
  mutate(pct = as.numeric(pct)) %>%
  replace_na(list("pct" = 0)) %>%
  mutate(normed_frac = pct/sum(pct)) 

pl_maj_summary <- pl_major %>% 
  #filter(strain %in% c("WT", "FM628")) %>%
  group_by(strain, pressure, class, replicate) %>%
  summarize(normed_frac = sum(normed_frac)) %>%
  group_by(strain, pressure, class) %>% 
  summarize(
    avg_frac = mean(normed_frac),
    sdev_frac = sd(normed_frac),
    serr_frac = sdev_frac/sqrt(n())
  )

library(dplyr)
library(tidyr)
library(broom)

# --- run t-test for each class and strain ---
ttest_results <- pl_major %>%
  filter(strain == "WT") %>%
  group_by(strain, class) %>%
  summarise(
    ttest = list(
      t.test(
        normed_frac[pressure == 1],
        normed_frac[pressure == 250],
        var.equal = TRUE        # or FALSE if variances differ
      ) %>% tidy()
    ),
    .groups = "drop"
  ) %>% 
  unnest(ttest) 

ttest_results


########################################################################################################################
###################################################### Figure 2A #######################################################

classes_to_plot <- c( "PE", "PC", "PA", "PS", "PI","PG", "LPE", "LPC", "LPI", "LPS", "LPA")

class_cols <- c(
  "PS" = "#E79AA0",
  "PE" = "#3AA655",
  "PC" = "#6FA8DC",
  "PI" = "#8E9AD6",
  "PG" = "#EDE684",
  "PA" = "#E5C15A"
)

# -----------------------------
# 1) Start from your parsed long table
#    Must contain: strain, pressure (optional), replicate, class, pct,
#                  chain1_length, chain2_length
# -----------------------------
df0 <- data_long %>%
  mutate(
    pct = as.numeric(pct),
    chain1_length = readr::parse_number(as.character(chain1_length)),
    chain2_length = readr::parse_number(as.character(chain2_length))
  )

# -----------------------------
# 2) Keep only classes you want to plot
# -----------------------------
df1 <- df0 %>%
  filter(class %in% classes_to_plot)

# -----------------------------
# 3) Filter out odd-chain PE species ONLY (keep even-even PE)
#    Non-PE classes are untouched.
# -----------------------------
df2 <- df1 %>%
  filter(
    class != "PE" |
      (
        !is.na(chain1_length) & !is.na(chain2_length) &
          chain1_length %% 2 == 0 &
          chain2_length %% 2 == 0
      )
  )

# -----------------------------
# 4) Sum species -> class within each sample
# -----------------------------
# Define a "sample unit" (edit if you want to include pressure in x-axis)
class_by_sample <- df2 %>%
  group_by(strain, pressure, replicate, class) %>%
  summarise(class_pct = sum(pct, na.rm = TRUE), .groups = "drop")

# -----------------------------
# 5) Normalize within each sample across ONLY the plotted classes
# -----------------------------
class_by_sample <- class_by_sample %>%
  group_by(strain, pressure, replicate) %>%
  mutate(class_frac = class_pct / sum(class_pct, na.rm = TRUE)) %>%
  ungroup()

# -----------------------------
# 6) Mean ± SEM across replicates (within strain; optionally within pressure too)
#    If your example is single-pressure, you can filter pressure first or average over it.
# -----------------------------
class_summary <- class_by_sample %>%
  group_by(strain, class) %>%   # add pressure here if you want separate pressure panels
  summarise(
    mean = mean(class_frac, na.rm = TRUE)*100,
    sem  = sd(class_frac, na.rm = TRUE) / sqrt(sum(!is.na(class_frac)))*100,
    .groups = "drop"
  ) %>%
  mutate(
    class  = factor(class, levels = classes_to_plot),
    strain = factor(strain, levels = unique(strain))
  )


class_summary <- class_summary %>%
  mutate(
    strain = factor(strain,
                    levels = c("W303", "pem2", "CHO", "psd1"))
  )

# -----------------------------
# 7) Plot like your example: one panel per class, free y, error bars
# -----------------------------
ggplot(class_summary, aes(x = strain, y = mean, fill = class)) +
  geom_col(width = 0.7, alpha = 0.9) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                width = 0.15, linewidth = 0.6) +
  facet_wrap(~class, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = class_cols) +
  labs(x = NULL, y = "Mole Fraction") +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(1.2, "lines")
  )



########################################################################################
#######################################  Figure 2C PE/PC growth ########################################################
###############################################################################################
# Load the data
data <- read.csv("~/PC_PE_growth.csv")

# Rename columns to match the correct format and pivot the data
data_long <- data %>%
  rename(`1 bar` = X1.bar, `250 bar` = X250.bar) %>%
  pivot_longer(cols = c("1 bar", "250 bar"), names_to = "Pressure", values_to = "OD") %>%
  filter(!is.na(strain))  # Keep this to remove NA strains

# Set the custom order for strains
data_long$strain <- factor(data_long$strain, levels = c("W303", "pem2", "CHO", "psd1"))

# Calculate mean and SD
data_summary <- data_long %>%
  group_by(strain, Pressure) %>%
  summarise(mean_OD = mean(OD, na.rm = TRUE), sd_OD = sd(OD, na.rm = TRUE), .groups = "drop")

#grouped bar plot
ggplot(data_summary, aes(x = strain, y = mean_OD, fill = Pressure)) +
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.75)) +
  geom_errorbar(
    aes(ymin = mean_OD - sd_OD, ymax = mean_OD + sd_OD),
    position = position_dodge(width = 0.75),
    width = 0.2,
    color = "black"
  ) +
  geom_jitter(
    data = data_long,
    aes(x = strain, y = OD, fill = Pressure),
    color = "black",
    size = 0.5,
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75)
  ) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  labs(x = NULL, y = expression(Delta~OD[600]), fill = "Pressure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5)) +
  theme(legend.position = 'none')


#########################################################################
#growth ratio
growth_ratio <- data %>%
  filter(!is.na(strain)) %>%
  mutate(ratio = X250.bar / X1.bar)

# Summarize by strain
growth_ratio_summary <- growth_ratio %>%
  group_by(strain) %>%
  summarise(
    mean_ratio = mean(ratio, na.rm = TRUE),
    sd_ratio = sd(ratio, na.rm = TRUE),
    se_ratio = sd_ratio / sqrt(n()),
    .groups = "drop"
  )

# Order strains for plotting
growth_ratio_summary$strain <- factor(growth_ratio_summary$strain,
                                      levels = c("W303", "pem2", "CHO", "psd1"))

ggplot(growth_ratio_summary, aes(x = strain, y = mean_ratio, fill = strain)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = mean_ratio - se_ratio, ymax = mean_ratio + se_ratio),
                width = 0.2, color = "black") +
  scale_fill_grey(start = 0.3, end = 0.8) +  # grayscale from darker to lighter
  labs(x = NULL, y = "Relative growth") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  scale_y_continuous(limits = c(0, 0.75), breaks = seq(0, 0.75, 0.25))

# choose your control
growth_ratio <- growth_ratio %>%
  mutate(strain = factor(strain)) %>%
  mutate(strain = relevel(strain, ref = "W303"))

# one-way ANOVA
fit <- aov(ratio ~ strain, data = growth_ratio)
summary(fit)

# Dunnett contrasts (each strain vs control)
dun <- glht(fit, linfct = mcp(strain = "Dunnett"))
sdun <- summary(dun)

# table of Dunnett results
dun_tab <- data.frame(
  contrast = names(sdun$test$coefficients),
  estimate = unname(sdun$test$coefficients),
  t       = unname(sdun$test$tstat),
  p.value = unname(sdun$test$pvalues),
  row.names = NULL
)

# make sure dun_tab has 'strain' and 'sig'
dun_tab <- dun_tab %>%
  mutate(strain = sub(" - .*", "", contrast),
         sig = cut(p.value,
                   breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                   labels = c("***","**","*",".","ns")))

# summarise strain means
sum_tab <- growth_ratio %>%
  group_by(strain) %>%
  summarise(mean = mean(ratio, na.rm = TRUE),
            sem  = sd(ratio, na.rm = TRUE)/sqrt(dplyr::n()),
            .groups = "drop")

# join with Dunnett results
plot_tab <- sum_tab %>%
  left_join(dplyr::select(dun_tab, strain, sig), by = "strain") %>%
  mutate(sig = ifelse(is.na(sig), "", as.character(sig)))


# plot with stars above points
ggplot(plot_tab, aes(x = strain, y = mean)) +
  geom_col(width = 0.7, fill = "grey80") +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
  geom_text(aes(y = mean + sem, label = sig), vjust = -0.5, size = 5) +
  labs(y = "Growth ratio (250/1 bar)", x = NULL,
       title = "Growth ratios by strain (Dunnett vs W303)") +
  theme_minimal(base_size = 12)



##########################################################################################################################
####################################################  Supp. Fig 1A   #######################################################
##########################################################################################################################  

# 1) PE chains as before, using normed_frac which is normalized to total PL pool
pe_chains <- pl_major %>%
  filter(class == "PE") %>%
  mutate(
    chain1_length = readr::parse_number(as.character(chain1_length)),
    chain2_length = readr::parse_number(as.character(chain2_length))
  ) %>%
  pivot_longer(
    cols = c(chain1_length, chain2_length),
    names_to = "position",
    values_to = "chain_length"
  ) %>%
  filter(!is.na(chain_length), chain_length > 0) %>%
  mutate(
    parity = ifelse(chain_length %% 2 == 0, "Even", "Odd"),
    chain_contrib = normed_frac / 2
  )

# 2) Sum odd/even contributions per sample
pe_parity_abs <- pe_chains %>%
  group_by(strain, pressure, replicate, parity) %>%
  summarise(total = sum(chain_contrib), .groups = "drop")
# NOTE: no renormalization step here

# 3) Mean ± SEM across replicates
pe_parity_abs_summary <- pe_parity_abs %>%
  group_by(strain, pressure, parity) %>%
  summarise(
    mean = mean(total, na.rm = TRUE),
    sem  = sd(total, na.rm = TRUE) / sqrt(sum(!is.na(total))),
    .groups = "drop"
  )

# 4) Plot: absolute mole fraction contribution (relative to total pool)
ggplot(pe_parity_abs_summary,
       aes(x = strain, y = mean, fill = parity)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, color = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                position = position_dodge(width = 0.7),
                width = 0, linewidth = 0.4) +
  facet_wrap(~pressure) +
  scale_fill_manual(values = c("Even" = "grey40", "Odd" = "#D55E00")) +
  labs(
    x = NULL,
    y = "PE (mole fraction)",
    fill = "Chain parity"
  ) +
  theme_classic(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################################################################################################
############################################# Figure 4B and Supp Fig 3I #################################

#column bars by class vs pressure
#YES

key_strain = c(
  "WT" = "W303",
  "FM628" = "FM628"
)
pl_maj_summary <- pl_maj_summary %>%
  mutate(class = factor(class, levels = c("PS", "PA", "LPC", "PC", "LPE", "PE", "LPI", "PI"  
                                          ))) 


# Dodged barplot with points and error bars
pl_maj_summary %>% 
  filter(strain == "WT") %>%
  ggplot() +
  facet_wrap(~class, nrow = 1, ncol = 11, scales = "free_y") +
  geom_col(
    #data = pl_maj_summary, 
    position = position_dodge(width = 225),
    aes(
      x = pressure,
      y = avg_frac*100,
      fill = class
    )
  ) + 
  geom_errorbar(
    #data = pl_maj_summary %>%  filter(strain == "WT"),
    position = position_dodge(width = 225), width = 0, 
    aes(
      x = pressure,
      ymin = (avg_frac - serr_frac)*100,
      ymax = (avg_frac + serr_frac)*100,
      group = class
    )
  ) +
  scale_x_continuous(breaks = c(1, 250)) +
  theme_pubr(base_size = 6) +
  theme(legend.position="none") +
  labs(
    x = "Pressure (bar)",
    y = "Mol % phospholipids"
  ) +
  scale_fill_manual(values = chroma_cl) 



pl_maj_summary %>% 
  filter(strain == "FM628") %>%
  ggplot(aes(x = factor(pressure))) +  # convert pressure to factor for dodging
  facet_wrap(~class, ncol=2, nrow=3, scales = "free_y") +
  geom_line(
    aes(y = avg_frac * 100, group = class, color = class),
    position = position_dodge(width = 0.5), 
    size=1
  ) +
  geom_bar(
    aes(
      ymin = (avg_frac - serr_frac) * 100,
      ymax = (avg_frac + serr_frac) * 100,
      group = class
    ),
    position = position_dodge(width = 0.5),
    width = 0.3, size = 0.5
  ) +
  geom_point(
    data = pl_major %>%
      filter(strain == "WT") %>%
      group_by(class, pressure, strain, replicate) %>%
      summarize(normed_frac = sum(normed_frac) * 100),
    aes(x = factor(pressure), y = normed_frac, group = class),
    position = position_dodge(width = 0.5), 
    size = 0.5
  ) +
  scale_x_discrete(name = "Pressure (bar)") +
  theme_pubr(base_size = 6) +
  theme(legend.position = "none") +
  labs(y = "% in Lipid Class")

##########################################################################
############################################# Figure 5A and Supp Fig 3H ######################## 
##########################################################################################################################################
#stacked bar plot of PL composition per strain and pressure

# Define the key_strain mapping
#YES
key_strain <- c("WT" = "W303", "FM628" = "FM628")

#key_strain <- c("W303" = "W303", "pem2" = "pem2", "psd1" = "psd1", "CHO" = "CHO")

# Modify the strain column to ensure the desired order of facets
pl_maj_summary <- pl_maj_summary %>% 
  filter(strain %in% c("WT", "FM628")) %>%
  #c("W303", "pem2", "psd1", "CHO")) %>%
  mutate(strain = factor(strain, levels = names(key_strain))) %>%
  filter(class %in% c("DAG"))

# Plot with the customized facet order
ggplot() +
  facet_wrap(~ strain, ncol = 5, labeller = labeller(strain = key_strain)) +
  labs(x = "Pressure (bar)", y = "mole fraction") +
  scale_x_continuous(breaks = pl_maj_summary %>% pull(pressure) %>% unique) +
  geom_col(
    data = pl_maj_summary, 
    aes(
      x = pressure,
      y = avg_frac,
      fill = class
      #fill = factor(class, levels = c("PS", "PG", "PC", "PA", "PI", "PE"))
    )
  ) +
  geom_errorbar(data = pl_maj_summary, 
                aes(
                  x = pressure,
                  ymin = avg_frac - serr_frac, 
                  ymax = avg_frac + serr_frac, 
                  group = strain),
                width = 0.15, linewidth = 0.6) +
  labs(fill = "Lipid") +
  theme_pubclean(base_size = 6) +
  #scale_fill_manual(values = chroma_cl) +
  #scale_fill_locuszoom() +
  theme(legend.position = "none") #+ 
theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())


########################################
# Supp Fig 3C, D
#stacked db in class vs pressure
#YES

key_strain = c(
  "WT" = "W303",
  "FM628" = "FM628"
)


# Filter and normalize data within each pressure condition and class
pl_major <- data_long %>%
  filter(class %in% c("PE", "PC", "PA", "PS", "PI", "DAG")) %>%
  mutate(pct = as.numeric(pct)) %>%
  replace_na(list(pct = 0))

# Summarize data
summary_data <- pl_major %>%
  filter(strain %in% c("FM628")) %>%
  #filter(total_DB <= 2) %>%
  group_by(strain, pressure, class, total_DB, replicate) %>%
  summarize(normed_frac = sum(pct), .groups = 'drop') %>%
  group_by(strain, pressure, class, replicate) %>%
  mutate(normed_frac = normed_frac / sum(normed_frac)) %>%
  group_by(strain, pressure, class, total_DB) %>%
  summarize(
    avg_frac = mean(normed_frac),
    sdev_frac = sd(normed_frac),
    serr_frac = sdev_frac / sqrt(n()),
    .groups = 'drop'
  )

dark2_palette <- brewer.pal(5, "Dark2")
names(dark2_palette) <- c(4, 3, 2, 1, 0)

# Plotting the data
ggplot(data = summary_data, aes(x = factor(pressure), y = avg_frac*100, fill = factor(total_DB))) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ class, scales = "free_y") +
  labs(
    x = "Pressure (bar)",
    y = "mole percent (%)",
    fill = "Total DB"
    #title = "Normalized Mole Percent of Lipid Classes by Pressure and Strain WT"
  ) +
  theme_minimal() +
  scale_fill_manual(values = dark2_palette) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )


#######################
#Supp Fig 3B
#total db plot with error bars final
#######################
db_maj_summary <- pl_major %>% 
  filter(strain %in% c("WT", "FM628"))  %>%
  filter(total_DB >= 0 & total_DB <= 4) %>%
  group_by(strain, pressure, total_DB, replicate) %>%
  summarize(db_frac = sum(pct)) %>%
  group_by(strain, pressure, total_DB) %>% 
  summarize(
    avg_frac = mean(db_frac),
    sdev_frac = sd(db_frac),
    serr_frac = sdev_frac/sqrt(n()),
  )

t_test_results <- pl_major %>%
  filter(total_DB >= 0 & total_DB <= 4) %>%
  group_by(strain, pressure, total_DB, replicate) %>%
  summarize(db_frac = sum(pct)) %>%
  group_by(strain, total_DB) %>%
  summarize(
    p_value = t.test(db_frac ~ pressure)$p.value
  )

dark2_palette <- brewer.pal(3, "Dark2")
names(dark2_palette) <- c(0, 1, 250)

db_maj_summary %>% 
  filter(total_DB %in% (0:4)) %>% ###double check filtering is with %in%
  ggplot() +
  facet_wrap(~factor(key_strain[strain], levels = c("W303", "FM628")), ncol = 3) + ## changed this 
  geom_col(
    #data = pl_maj_summary, 
    position = position_dodge(),
    aes(
      x = total_DB,
      y = avg_frac*100,
      fill = as.factor(pressure)
    )
  ) + 
  #facet_grid(cols = vars(key_strain[strain])) + ###commented whis out
  labs(
    x = "Total Double Bonds",
    y = "mol % phospolipids", 
    fill = "Pressure (bar)"
  ) +
  #theme_classic2()+
  scale_fill_manual(values = dark2_palette) +
  theme(legend.position = "right")+ 
  geom_errorbar(
    position = position_dodge(), 
    aes(
      x = total_DB,
      ymin = (avg_frac - serr_frac)*100,
      ymax = (avg_frac + serr_frac)*100,
      group = pressure
    )
  ) +
  geom_point(
    data = pl_major %>%
      filter(strain %in% c("FM628", "WT"))  %>%
      filter(total_DB %in% (0:4)) %>%
      group_by(total_DB, strain, replicate, pressure) %>%
      summarize(normed_frac = sum(normed_frac) * 100),
    position = position_dodge(width = 0.9),
    aes(
      x = total_DB,
      y = normed_frac,
      group = pressure
    )
  ) +
  theme_minimal() +
  theme(legend.position = "none")



#########################################
#sn1 db
#YES
##############################################

db_maj_summary <- pl_major %>% 
  filter(strain %in% c("WT", "FM628"))  %>%
  filter(class %in% c('PC', 'PE', "PA", 'PI', "PS", "PG")) %>%
  filter(total_DB >= 0 & total_DB <= 4) %>%
  group_by(strain, pressure, chain1_db, replicate) %>%
  summarize(db_frac = sum(normed_frac)) %>%
  group_by(strain, pressure, chain1_db) %>% 
  summarize(
    avg_frac = mean(db_frac),
    sdev_frac = sd(db_frac),
    serr_frac = sdev_frac/sqrt(n()),
  )

t_test_results <- pl_major %>%
  filter(total_DB >= 0 & total_DB <= 4) %>%
  group_by(strain, pressure, chain1_db, replicate) %>%
  summarize(db_frac = sum(normed_frac)) %>%
  group_by(strain, chain1_db) %>%
  summarize(
    p_value = t.test(db_frac ~ pressure)$p.value
  )

dark2_palette <- brewer.pal(3, "Dark2")
names(dark2_palette) <- c(0, 1, 250)

db_maj_summary %>% 
  filter(chain1_db %in% (0:4)) %>% ###double check filtering is with %in%
  ggplot() +
  facet_wrap(~factor(key_strain[strain], levels = c("W303", "FM628")), ncol = 3) + ## changed this 
  geom_col(
    #data = pl_maj_summary, 
    position = position_dodge(),
    aes(
      x = chain1_db,
      y = avg_frac*100,
      fill = as.factor(pressure)
    )
  ) + 
  #facet_grid(cols = vars(key_strain[strain])) + ###commented whis out
  labs(
    x = "SN1 Chain Double Bonds",
    y = "mol % phospolipids", 
    fill = "Pressure (bar)"
  ) +
  #theme_classic2()+
  scale_fill_manual(values = dark2_palette) +
  theme(legend.position = "right")+ 
  geom_errorbar(
    position = position_dodge(), 
    aes(
      x = chain1_db,
      ymin = (avg_frac - serr_frac)*100,
      ymax = (avg_frac + serr_frac)*100,
      group = pressure
    )
  ) +
  geom_point(
    data = pl_major %>%
      filter(strain %in% c("FM628", "WT"))  %>%
      filter(chain1_db %in% (0:4)) %>%
      group_by(chain1_db, strain, replicate, pressure) %>%
      summarize(normed_frac = sum(normed_frac) * 100),
    position = position_dodge(width = 0.9),
    aes(
      x = chain1_db,
      y = normed_frac,
      group = pressure
    )
  ) +
  theme(legend.position = "none")


##################################################
#sn2 db
#YES 

db_maj_summary <- pl_major %>% 
  filter(strain %in% c("WT", "FM628"))  %>%
  filter(class %in% c('PC', 'PE', "PA", 'PI', "PS")) %>%
  filter(total_DB >= 0 & total_DB <= 4) %>%
  group_by(strain, pressure, chain2_db, replicate) %>%
  summarize(db_frac = sum(normed_frac)) %>%
  group_by(strain, pressure, chain2_db) %>% 
  summarize(
    avg_frac = mean(db_frac),
    sdev_frac = sd(db_frac),
    serr_frac = sdev_frac/sqrt(n()),
  )

t_test_results <- pl_major %>%
  filter(total_DB >= 0 & total_DB <= 4) %>%
  group_by(strain, pressure, chain2_db, replicate) %>%
  summarize(db_frac = sum(normed_frac)) %>%
  group_by(strain, chain2_db) %>%
  summarize(
    p_value = t.test(db_frac ~ pressure)$p.value
  )

dark2_palette <- brewer.pal(3, "Dark2")
names(dark2_palette) <- c(0, 1, 250)

db_maj_summary %>% 
  filter(chain2_db %in% (0:3)) %>% ###double check filtering is with %in%
  ggplot() +
  facet_wrap(~factor(key_strain[strain], levels = c("W303", "FM628")), ncol = 3) + ## changed this 
  geom_col(
    #data = pl_maj_summary, 
    position = position_dodge(),
    aes(
      x = chain2_db,
      y = avg_frac*100,
      fill = as.factor(pressure)
    )
  ) + 
  #facet_grid(cols = vars(key_strain[strain])) + ###commented whis out
  labs(
    x = "SN2 chain double bonds",
    y = "mol % phospolipids", 
    fill = "Pressure (bar)"
  ) +
  #theme_classic2()+
  scale_fill_manual(values = dark2_palette) +
  theme(legend.position = "right")+ 
  geom_errorbar(
    position = position_dodge(), 
    aes(
      x = chain2_db,
      ymin = (avg_frac - serr_frac)*100,
      ymax = (avg_frac + serr_frac)*100,
      group = pressure
    )
  ) +
  geom_point(
    data = pl_major %>%
      filter(strain %in% c("FM628", "WT"))  %>%
      filter(chain2_db %in% (0:3)) %>%
      group_by(chain2_db, strain, replicate, pressure) %>%
      summarize(normed_frac = sum(normed_frac) * 100),
    position = position_dodge(width = 0.9),
    aes(
      x = chain2_db,
      y = normed_frac,
      group = pressure
    )
  ) +
  theme(legend.position = "none")


#############################
# Supp Fig 3B
########total length###########
###############################

len_maj_summary <- pl_major %>% 
  filter(strain %in% c("WT", "FM628"))  %>%
  filter(total_length >= 29 & total_length <= 36) %>%
  group_by(strain, pressure, total_length, replicate) %>%
  summarize(db_frac = sum(pct)) %>%
  group_by(strain, pressure, total_length) %>% 
  summarize(
    avg_frac = mean(db_frac),
    sdev_frac = sd(db_frac),
    serr_frac = sdev_frac/sqrt(n())
  )

t_test_results <- pl_major %>%
  filter(total_length >= 29 & total_length <= 36) %>%
  group_by(strain, pressure, total_length, replicate) %>%
  summarize(len_frac = sum(pct)) %>%
  group_by(strain, total_length) %>%
  summarize(
    p_value = t.test(len_frac ~ pressure)$p.value
  )

dark2_palette <- brewer.pal(3, "Dark2")
names(dark2_palette) <- c(0, 1, 250)


len_maj_summary %>% 
  ggplot() +
  #facet_wrap(~strain, ncol=3) +
  facet_wrap(~factor(key_strain[strain], levels = c("W303", "FM628")), ncol = 3) + ## changed this 
  geom_col(
    #data = pl_maj_summary, 
    position = position_dodge(),
    aes(
      x = total_length,
      y = avg_frac*100,
      fill = as.factor(pressure)
    )
  ) + 
  #facet_grid(cols = vars(key_strain[strain])) +
  labs(
    x = "Total length",
    y = "mol % phospholipids", 
    fill = "Pressure (bar)"
  ) +   
  scale_fill_manual(values = dark2_palette) +
  #theme_classic2()+
  theme(legend.position = "none")+ 
  geom_errorbar(
    position = position_dodge(), 
    aes(
      x = total_length,
      ymin = (avg_frac - serr_frac)*100,
      ymax = (avg_frac + serr_frac)*100,
      group = pressure
    )
  ) +
  geom_point(
    data = pl_major %>%
      filter(strain %in% c("FM628", "WT"))  %>%
      filter(total_length >= 29 & total_length <= 36) %>%
      group_by(total_length, strain, replicate, pressure) %>%
      summarize(normed_frac = sum(normed_frac) * 100),
    position = position_dodge(width = 0.9),
    aes(
      x = total_length,
      y = normed_frac,
      group = pressure
    )
  )


#############################################
############################################
#total length by class
#Supp Fig 3E-F
# Filter and normalize data within each pressure condition and class
pl_major <- data_long %>%
  filter(class %in% c("PE", "PC", "PA", "PS", "PI", "DAG")) %>%
  mutate(pct = as.numeric(pct)) %>%
  replace_na(list(pct = 0))

# Summarize data
summary_data <- pl_major %>%
  filter(strain == "FM628") %>%
  filter(total_length %% 2 == 0) %>%
  group_by(strain, pressure, class, total_length, replicate) %>%
  summarize(normed_frac = sum(pct), .groups = 'drop') %>%
  group_by(strain, pressure, class, replicate) %>%
  mutate(normed_frac = normed_frac / sum(normed_frac)) %>%
  group_by(strain, pressure, class, total_length) %>%
  summarize(
    avg_frac = mean(normed_frac),
    sdev_frac = sd(normed_frac),
    serr_frac = sdev_frac / sqrt(n()),
    .groups = 'drop'
  )

dark2_palette <- brewer.pal(5, "Dark2")
names(dark2_palette) <- c(4, 3, 2, 1, 0)

# Plotting the data
ggplot(data = summary_data, aes(x = factor(pressure), y = avg_frac*100, fill = factor(total_length))) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ class, scales = "free_y") +
  labs(
    x = "Pressure (bar)",
    y = "mole % lipid",
    fill = "Total length"
    #title = "Normalized Mole Percent of Lipid Classes by Pressure and Strain WT"
  ) +
  theme_minimal() +
  #scale_fill_manual(values = dark2_palette) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

########chain 1 comp##################
#######################################

key_strain = c(
  "WT" = "W303",
  "FM628" = "FM628"
)

ch1db_maj_summary <- pl_major %>% 
  filter(chain2_db >= 0 & chain2_db <= 4) %>%
  filter(chain1_length %in% c(16,18)) %>%
  group_by(strain, pressure, chain1, replicate) %>%
  summarize(db_frac = sum(normed_frac)) %>%
  group_by(strain, pressure, chain1) %>% 
  summarize(
    avg_frac = mean(db_frac),
    sdev_frac = sd(db_frac),
    serr_frac = sdev_frac/sqrt(n())
  ) 

t_test_results <- pl_major %>%
  filter(chain2_db >= 0 & chain2_db <= 4) %>%
  filter(chain1_length %in% c(16,18)) %>%
  group_by(strain, pressure, chain1, replicate) %>%
  summarize(db_frac = sum(normed_frac)) %>%
  group_by(strain, chain1) %>%
  summarize(
    p_value = t.test(db_frac ~ pressure)$p.value
  )

dark2_palette <- brewer.pal(3, "Dark2")
names(dark2_palette) <- c(0, 1, 250)

ch1db_maj_summary %>% 
  #filter(avg_frac >= 0.0025) %>%
  ggplot() +
  #facet_wrap(~strain, ncol=3) +
  facet_wrap(~factor(key_strain[strain], levels = c("W303", "FM628")), ncol = 3) + ## changed this 
  geom_col(
    #data = pl_maj_summary, 
    position = position_dodge(),
    aes(
      x = chain1,
      y = avg_frac*100,
      fill = as.factor(pressure)
    )
  ) + 
  #facet_grid(cols = vars(key_strain[strain])) +
  labs(
    x = "Chain sn1",
    y = "mol % phospholipids", 
    fill = "Pressure (bar)"
  ) +
  #theme_classic2()+
  theme(legend.position = "none")+ 
  geom_errorbar(
    position = position_dodge(), 
    aes(
      x = chain1,
      ymin = (avg_frac - serr_frac)*100,
      ymax = (avg_frac + serr_frac)*100,
      group = pressure
    )
  ) +
  scale_fill_manual(values = dark2_palette) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_point(
    data = pl_major %>%
      filter(chain1_db >= 0 & chain2_db <= 3) %>%
      filter(chain1_length %in% c(16,18)) %>%
      group_by(chain1, strain, replicate, pressure) %>%
      summarize(normed_frac = sum(normed_frac) * 100),
    position = position_dodge(width = 0.9),
    aes(
      x = chain1,
      y = normed_frac,
      group = pressure
    )
  ) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")


########################################
########chain 2 db##################
#######################################

key_strain = c(
  "FM628" = "FM628",
  "WT" = "W303"
)

ch2db_maj_summary <- pl_major %>% 
  filter(chain2_db >= 0 & chain2_db <= 3) %>%
  filter(chain2_length %in% c(16,18)) %>%
  group_by(strain, pressure, chain2, replicate) %>%
  summarize(db_frac = sum(normed_frac)) %>%
  group_by(strain, pressure, chain2) %>% 
  summarize(
    avg_frac = mean(db_frac),
    sdev_frac = sd(db_frac),
    serr_frac = sdev_frac/sqrt(n())
  ) 

t_test_results <- pl_major %>%
  filter(chain2_db >= 0 & chain2_db <= 4) %>%
  filter(chain2_length %in% c(16,18)) %>%
  group_by(strain, pressure, chain2, replicate) %>%
  summarize(db_frac = sum(normed_frac)) %>%
  group_by(strain, chain2) %>%
  summarize(
    p_value = t.test(db_frac ~ pressure)$p.value
  )

dark2_palette <- brewer.pal(3, "Dark2")
names(dark2_palette) <- c(0, 1, 250)

ch2db_maj_summary %>% 
  #filter(strain %in% c("WT"))  %>%
  #filter(avg_frac >= 0.0025) %>%
  ggplot() +
  #facet_wrap(~strain, ncol=3) +
  facet_wrap(~factor(key_strain[strain], levels = c("W303", "FM628")), ncol = 3) + ## changed this 
  geom_col(
    #data = ch2db_maj_summary, 
    position = position_dodge(),
    aes(
      x = chain2,
      y = avg_frac*100,
      fill = as.factor(pressure)
    )
  ) + 
  #facet_grid(cols = vars(key_strain[strain])) +
  labs(
    x = "Chain sn2",
    y = "mol % phospholipids", 
    fill = "Pressure (bar)"
  ) +
  #theme_classic2()+
  geom_errorbar(
    position = position_dodge(), 
    aes(
      x = chain2,
      ymin = (avg_frac - serr_frac)*100,
      ymax = (avg_frac + serr_frac)*100,
      group = pressure
    )
  ) +
  scale_fill_manual(values = dark2_palette) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_point(
    data = pl_major %>%
      filter(chain2_db >= 0 & chain2_db <= 3) %>%
      filter(chain2_length %in% c(16,18)) %>%
      group_by(chain2, strain, replicate, pressure) %>%
      summarize(normed_frac = sum(normed_frac) * 100),
    position = position_dodge(width = 0.9),
    aes(
      x = chain2,
      y = normed_frac,
      group = pressure
    )
  ) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")
  


##########################################################################
###########################################################################
###########################################################################
##### 2chain composition vs pressure
#YES

key_strain = c(
  "FM628" = "FM628"
  "WT" = "W303"
)

chains_maj_summary <- pl_major %>% 
  filter(chain2_db <= 3) %>%
  filter(chain2_length < 20) %>%
  filter(chain2_length %% 2 == 0) %>%
  filter(total_length %% 2 == 0) %>%
  group_by(strain, pressure, chain, replicate) %>%
  summarize(db_frac = sum(normed_frac)) %>%
  group_by(strain, pressure, chain) %>% 
  summarize(
    avg_frac = mean(db_frac),
    sdev_frac = sd(db_frac),
    serr_frac = sdev_frac/sqrt(n()), 
    .groups = 'drop'
  ) 


dark2_palette <- brewer.pal(3, "Dark2")
names(dark2_palette) <- c(0, 1, 250)

chains_maj_summary %>% 
  filter(strain %in% c("WT"))  %>%
  ggplot() +
  #facet_wrap(~strain, ncol=3) +
  #facet_wrap(~factor(key_strain[strain], levels = c("W303", "FM628")), ncol = 2) + ## changed this 
  geom_col(
    #data = ch2db_maj_summary, 
    position = position_dodge(),
    aes(
      x = chain,
      y = avg_frac*100,
      fill = as.factor(pressure)
    )
  ) + 
  #facet_grid(cols = vars(key_strain[strain])) +
  labs(
    x = "Acyl chains",
    y = "mol % phospholipids", 
    fill = "Pressure (bar)"
  ) +
  #theme_classic2()+
  theme(legend.position = "none")+ 
  geom_errorbar(
    position = position_dodge(), 
    aes(
      x = chain,
      ymin = (avg_frac - serr_frac)*100,
      ymax = (avg_frac + serr_frac)*100,
      group = pressure
    )
  ) +
  scale_fill_manual(values = dark2_palette) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) #+
  geom_point(
    data = pl_major %>%
      filter(chain2_db >= 0 & chain2_db <= 3) %>%
      #filter(chain2_length %in% c(16,18)) %>%
      group_by(chain, strain, replicate, pressure) %>%
      summarize(normed_frac = sum(normed_frac) * 100),
    position = position_dodge(width = 0.9),
    aes(
      x = chain,
      y = normed_frac,
      group = pressure
    )
  ) 
  



################################################################################################
##################################################################################################
##### FM628 vs W303 fitness ##########################
#Supp Fig 3A####################
################################################################################################
##################################################################################################

data <- read.csv("~/FM628_W303_growth.csv") 

# Rename columns to match the correct format and pivot the data
data_long <- data %>%
  rename (`1 bar` = X1.bar, `250 bar` = X250.bar) %>%  # Rename columns directly
  pivot_longer(cols = c("1 bar", "250 bar"), names_to = "Pressure", values_to = "OD") %>%
  filter(!is.na(strain)) %>%  # Filter out NA strains
  filter(Pressure == "250 bar") #filter which pressure condition to select. Can comment out for both together

# Set the custom order for strains
data_long$strain <- factor(data_long$strain, levels = c("W303", "FM628"))

# Calculate mean and SD
data_summary <- data_long %>%
  group_by(strain, Pressure) %>%
  summarise(mean_OD = mean(OD, na.rm = TRUE), sd_OD = sd(OD, na.rm = TRUE), .groups = "drop")

# Plot with facet grid in columns
ggplot(data_summary, aes(x = strain, y = mean_OD, fill = strain)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = mean_OD - sd_OD, ymax = mean_OD + sd_OD), width = 0.25, color = "black") +
  geom_jitter(data = data_long, aes(x = strain, y = OD), 
              color = "black", width = 0.25, height = 0, size = 1) + 
  scale_fill_manual(values = c("W303" = "gray", 
                               "FM628" = "gray35")) + 
  labs(x = "Strain", y = expression(Delta~OD600)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  facet_wrap(~ Pressure, nrow = 2) +  # Display facets in columns
  scale_y_continuous(limits = c(0, 0.75), breaks = seq(0, 0.75, 0.25)) 


library(tidyverse)
df <- read_csv("FM628_W303_growth.csv")

df_long <- df %>%
  pivot_longer(cols = c(`1bar`, `250bar`),
               names_to = "pressure",
               values_to = "growth")

df_long$strain   <- factor(df_long$strain, levels = c("W303", "FM628"))
df_long$pressure <- factor(df_long$pressure, levels = c("1bar", "250bar"))

summary_df <- df_long %>%
  group_by(strain, pressure) %>%
  summarise(
    mean = mean(growth),
    sem  = sd(growth)/sqrt(n()),
    .groups = "drop"
  )

ggplot(summary_df, aes(x = pressure, y = mean, fill = strain)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.9) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                position = position_dodge(width = 0.7),
                width = 0.15, size = 0.9) +
  geom_jitter(data = df_long,
              aes(x = pressure, y = growth, color = strain),
              position = position_jitterdodge(jitter.width = 0.08,
                                              dodge.width = 0.7),
              size = 1, alpha = 0.85,
              inherit.aes = FALSE) +
  scale_fill_manual(values = c("W303" = "grey30",
                               "FM628" = "grey70")) +
  scale_color_manual(values = c("W303" = "black",
                                "FM628" = "black")) +
  scale_y_continuous(
    limits = c(0, 1.2),
    breaks = seq(0, 1, by = 0.25),
    expand = c(0, 0)
  ) +
  labs(y = "Growth", x = "Pressure (bar)", fill = "Strain", color = "Strain") +
  theme_classic(base_size = 6)


################################################################################################