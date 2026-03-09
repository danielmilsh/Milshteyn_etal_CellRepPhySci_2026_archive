library(here)
library(tidyverse)
library(tidymodels)
library(httr)
library(ggplot2)
library(ggpubr)

source("~/saxs_helpers.R")


files_profs = list.files(here("05-yeast", "pressure"), pattern = "*.dat", full.names = TRUE)
file_metdata = here("metadata.csv")

# factor order for guest lipids
lipids_ord = c("1 bar", "250 bar")

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
  # ditch the degraded SOPPE sample
  #filter(ts == .12) %>%
  #filter(temp == 85) %>%
  #filter(press == 0) %>%
  filter(pdir == "up") %>%
  #filter(!(samp %in% paste("DM", str_pad(seq(11, 15), 4, pad='0'), sep=''))) %>% 
  # crop q-domain for normalization
  filter(between(q, 0.05, 0.35)) %>% 
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
lipid_order <- c("1 bar","250 bar")

##################################################
#####stop here for phase or continue for d space and c0
###################################################
## Fit HII characteristic spacing
# expected peak spacing ratio
# first element 1 corresponds to d
ratio = c(
  1,
  sqrt(3),
  2,
  sqrt(7),
  3,
  sqrt(11)
)

## test every q value as d 
firstorder = saxsmean %>% 
  # create interpolator function for each profile
  group_by(samp, lipid_guest, frac_guest, temp, press, ts, pdir) %>% 
  summarize(
    iqfun  = approxfun(q, iq) %>% list(),
    errfun = approxfun(q, err) %>% list()
  ) %>% 
  ## join them back to the data
  full_join(saxsmean, by = c("samp", "lipid_guest", "press", "temp", "ts", "pdir")) %>% 
  # calculate total amplitude of peaks using each q-value as proposed d
  rowwise() %>% 
  # use sum(q*iq) to avoid overweighting low-angle peaks
  mutate(peaksum = sum(q * iqfun(q * ratio), na.rm = TRUE)) %>% 
  group_by(lipid_guest, temp, press, ts, pdir) %>% 
  filter(peaksum == max(peaksum))

# calculate the peak coords for visual goodness-of-fit assessment
peaks = firstorder %>% 
  crossing(tibble(ratio = ratio)) %>% 
  group_by(samp, lipid_guest, frac_guest.x, temp, press, ts, pdir) %>% 
  mutate(
    q = q*ratio,
    iq = iqfun[[1]](q),
    err = errfun[[1]](q)
  ) %>% 
  # order the factor for plotting
  mutate(
    lipid_guest = lipid_guest %>% 
      factor(levels = lipids_ord)
  ) %>% 
  drop_na(lipid_guest)
#select(q, iq, err, ratio)

# plot the fits
saxsmean %>%
  filter(
    (press < 21) & 
    (temp == 85) &
    #!(temp %% 5) & # remove "odd" pressures
      (pdir == "up") &
      (ts == 0)
  ) %>% 
  #filter(lipid_guest != "SOPPE") %>% # filter out bad data for now
  ggplot(
    aes(
      x = q,
      y = iq
    )
  ) +
  facet_grid(
    rows = vars(lipid_guest), 
    cols = vars(press)
  ) +
  scale_y_log10() +
  scale_color_brewer(palette = "Set1") +
  geom_line(
    aes(
      group = samp
    )
  ) +
  geom_point(
    data = peaks %>% 
      filter(
        (pdir == "up") &
          (temp == 85) &
          (ratio < 2.5) &
          (press < 21)# remove "odd" pressures
      ),
    shape = 6,
    aes(color = lipid_guest),
  ) +
  labs(
    x = "q (1/Å)",
    y = "I(q) (arb. normalized)"
  ) +
  theme_pubr()
ggsave(here("pressureYeast_Peakspicked_05202025.pdf"), width = 3, height = 4)



spacings = peaks %>% 
  filter(ratio == 1) %>% 
  #filter(!(temp %% 5)) %>%
  filter(press < 41) %>%
  filter(temp == 85) %>%
  #left_join(DOPEcnpl, by = "press") %>%
  # calc d-spacing in Å and curvature c in Å^-1
  mutate(
    d = 2*pi/q,
    # lipid mass fraction, fully hydrated
    w = 0.754, # DOPE at 85°C per Tate and Gruner 1989, Table 2
    # calc Rw
    r = sqrt(2*d^2 * (1-w) / (pi*sqrt(3))),
    rnpl = r + 3.94,
    c = -1/r,
    cnpl = -1/rnpl 
  ) %>%
  ungroup() %>%
  select(lipid_guest, press, cnpl)

spacings[] <- lapply(spacings, function(x) {
  if (is.list(x)) {
    # Convert list elements to strings
    sapply(x, function(y) {
      if (is.function(y)) {
        return(NA)  # Replace functions with NA
      } else {
        return(toString(y))  # Convert non-function elements to string
      }
    })
  } else {
    x  # Keep non-list columns as is
  }
})

write.csv(spacings,"~/c0_pressureYeast.csv", row.names = FALSE)

##############################################
#Fig 3B - curvature of pressure yeast Hii
library(dplyr)
library(ggplot2)
library(broom)

# Fit separate GLMs for each lipid_guest
glm_preds <- spacings %>%
  group_by(lipid_guest) %>%
  do({
    model <- glm(cnpl ~ press, data = .)
    preds <- predict(model, se.fit = TRUE)
    
    tibble(
      press = .$press,
      cnpl_pred = preds$fit,
      se_cnpl_pred = preds$se.fit
    )
  }) %>%
  ungroup()

glm_preds_filtered <- glm_preds %>%
  filter(press %in% c(0, 20)) %>%
  mutate(
    lipid_guest = factor(lipid_guest, levels = c("250 bar", "1 bar"))
  )
ggplot(glm_preds_filtered, aes(x = factor(press), y = cnpl_pred)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = cnpl_pred - se_cnpl_pred,
      ymax = cnpl_pred + se_cnpl_pred
    ),
    width = 0.5
  ) +
  facet_wrap(~lipid_guest) +
  labs(
    x = "Pressure (bar)",
    y = expression(italic(c))
  ) +
  scale_y_reverse(
    limits = c(-0.028, -0.0405),
    breaks = c(-0.03, -0.04)
  ) +
  theme_pubr()

ggsave(here("pressureYeast_c0_20250610.pdf"), width = 6, height = 6, units = "cm", device = cairo_pdf)


