library(here)
library(tidyverse)
library(tidymodels)
library(httr)
library(ggplot2)
library(ggpubr)

source("~/saxs_helpers2.R")

files_profs = list.files(here("02-profiles"), pattern = "*.dat", full.names = TRUE)
file_metdata = here("metadata.csv")

# factor order for guest lipids
lipids_ord = c("W303a", "pem2-", "CHO", "psd1-")


# and load the data in
metadata = file_metdata %>% read_csv()

# read in SAXS data
# this step takes a hot minute,
# can be sped up by omitting files that won't be used
saxsdata = files_profs %>% 
  read_saxsdats() %>% 
  # parse/massage the fname parts
  mutate(
    frm = as.numeric(frm),
    rep = max(rep, frm)
  ) %>% 
  select(-typ, -ser, -frm) #%>% 
# remove trailing letter from samp if it is there
#mutate(samp = str_remove(samp, regex("(A|B)")))

# bind to metadata
saxsmeta = saxsdata %>% 
  left_join(metadata, by = "samp") %>% 
  group_by(samp, lipid_guest, frac_guest, ts, temp, press, q)

# normalize shots to the same amplitude
miniq = 0.1
maxiq = 1
saxsnorm = saxsmeta %>% 
  # can filter out data here
  filter(ts == 0) %>%
  #filter(pdir == "up") %>%
  # crop q-domain for normalization
  filter(between(q, 0.03, 0.35)) %>% 
  # normalize I(q)
  group_by(fname) %>% 
  mutate(
    iq = iq - min(iq),
    iq = iq/max(iq),
    iq = iq + miniq
  )

# average replicate shots
saxsmean = saxsnorm %>% 
  #filter(press == 0) %>% # there's some 60C in there!
  filter(pdir == "up") %>%
  mutate(rep = as.numeric(rep)) %>% 
  group_by(samp, lipid_guest, frac_guest, ts, temp, press, pdir, q) %>% 
  # take only the last 3 shots at each state
  filter(rep > (max(rep) - 3)) %>% 
  summarize(
    iq = mean(iq),
    err = sqrt(sum(err**2))/sqrt(n()),
    fname = list(fname)
  ) %>% 
  # order the factor for plotting
  mutate(
    lipid_guest = lipid_guest %>% 
      factor(levels = lipids_ord)
  )

# Define the order and colors for lipid_guest
lipid_order <- c("W303a", "pem2-", "CHO", "psd1-")


saxsmean %>%
  filter(
    #(temp %in% c(70, 75, 80, 85)) &
    #(press == 0)  &
    (temp == 85) &
      #(press <= 200) &
      #!(press %% 100) & # remove "odd" pressures
      (pdir == "up") 
    #(samp %in% c("D3M0002"))
  ) %>% 
  filter(lipid_guest != "NA") %>% # filter out bad data for now
  ggplot(
    aes(
      x = q,
      y = iq
    )
  ) +
  facet_grid(
    rows = vars(press*10),
    cols = vars(factor(lipid_guest, levels = lipid_order))
  ) +
  scale_y_log10() +
  scale_color_manual(values = lipid_colors) +
  geom_line(
    aes(
      group = samp,
      color = "#000000"
    )
  ) +
  labs(
    color = "Guest lipid",
    x = "q (1/Å)",
    y = "I(q) (arb. normalized)") +
  theme(
    #legend.position = c(0.3, 0.75),
    legend.position = "right",
    legend.background = element_blank()
  ) +
  theme_pubr(x.text.angle = 45) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(here("PEPC_Psweep_3500bar.pdf"), width = 6, height = 8)