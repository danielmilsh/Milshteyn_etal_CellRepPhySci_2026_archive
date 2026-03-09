library(here)
library(tidyverse)
library(tidymodels)
library(httr)
library(ggplot2)
library(ggpubr)

source("~/saxs_helpers.R")

#SAXS PI profiles and c0
#Fig 4C-E
#######################################################################################
################################################################################
####### Pt1 is phase, pt2 is dspace + c0 ############################
###############################################################
files_profs = list.files(here("PI-profiles"), pattern = "*.dat", full.names = TRUE) #PI pressure data
file_metdata = here("metadata.csv")

# factor order for guest lipids
lipids_ord = c("DOPE", "DOPI (10%)", "Soy PI (10%)", "DOPC (20%)", "POPC (10%)")

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
  filter(ts == 0.12) %>%
  #filter(temp == 35) %>%
  #filter(press == 0) %>%
  filter(pdir == "up") %>%
  #filter(!(samp %in% paste("DM", str_pad(seq(11, 15), 4, pad='0'), sep=''))) %>% 
  # crop q-domain for normalization
  filter(between(q, 0.05, 0.35)) %>% 
  # normalize I(q)
  group_by(samp) %>% 
  mutate(
    iq = iq - min(iq),
    iq = iq/max(iq),
    iq = iq + miniq
  )

## PEEK NORMALIZED PROFILES
saxsnorm %>% 
  #filter(samp %in% c("DM0013", "DM0006")) %>%
  #filter((press == 0)) %>% # simplifying for now
  #filter(temp == 85) %>%
  #filter(ts == 0) %>% 
  filter(pdir == "up") %>%
  # take just the last profile per condition
  group_by(samp, lipid_guest, temp, press, q) %>% 
  summarize(iq = last(iq)) %>% 
  # draw plot
  ggplot(
    aes(
      x = q,
      y = iq,
      color = lipid_guest
    )
  ) +
  geom_line() +
  scale_y_log10() +
  theme_pubr() +
  labs(
    x = "q (1/Å)",
    y = "I(q) (arb.)"
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


saxsmean %>%
  filter(
    (temp == 35) &
    #(temp %in% c(65, 75, 85)) &
    (press < 100) &
    !(press %% 20) & # remove "odd" pressures
    (pdir == "up")) %>%
  #(lipid_guest %in% c("1 bar", "250 bar"))) %>%
  filter(lipid_guest != "NA") %>% # filter out bad data for now
  ggplot(
    aes(
      x = q,
      y = iq#, 
      #fill = "#000000"
    )
  ) +
  facet_grid(
    rows = vars(lipid_guest), 
    cols = vars(press)) +
  scale_y_log10() +
  #scale_color_manual(values = lipid_colors) +
  geom_line(
    aes(
      group = samp
    )
  ) +
  labs(
    x = "q (1/Å)",
    y = "I(q) (arb. normalized)") +
  xlim(0.05, 0.25) +
  theme(
    #legend.position = c(0.3, 0.75),
    legend.position = "right",
    legend.background = element_blank()
  ) +
  theme_pubr(x.text.angle = 45) +
  theme_classic() +
  theme(legend.position = "none")

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
    #(press == 0) & 
    #(press <= 200) &
    !(press %% 20) & # remove "odd" pressures
      (pdir == "up") #&
    #(ts == 0)
  ) %>% 
  ggplot(
    aes(
      x = q,
      y = iq
    )
  ) +
  facet_grid(
    rows = vars(ts),
    cols = vars(press)
  ) +
  scale_y_log10() +
  scale_color_brewer(palette = "Set1") +
  geom_line(
    aes(
      group = samp,
      color = lipid_guest
    )
  ) +
  geom_point(
    data = peaks %>% 
      filter(
        (pdir == "up") &
          (press <= 200) &
          !(press %% 20)# remove "odd" pressures
      ),
    aes(color = lipid_guest),
    shape = 6#,
    #color = "black"
  ) +
  labs(
    title = "Profiles by pressure and [tricosene], upsweep only, max sum(q*peaks)",
    color = "Guest lipid (20 mol %)",
    x = "q (1/Å)",
    y = "I(q) (arb. normalized)"
  ) +
  theme_pubr()


## plot tubule R vs P for all compositions, 12% TS
# this ignores peak-pick errors due to lamellar phase transition.
spacings = peaks %>% 
  filter(pdir == "up") %>%
  filter(ratio == 1) %>% 
  filter (temp == 35) %>%
  filter(!(press %% 20)) %>%
  #filter(press == 0) %>%
  # calc d-spacing in Å and curvature c in Å^-1
  mutate(
    d = 2*pi/q,
    # lipid mass fraction, fully hydrated
    w = 0.718, # DOPE at 35°C per Tate and Gruner 1989, Table 2
    # calc Rw
    r = sqrt(2*d^2 * (1-w) / (pi*sqrt(3))),
    rnpl = r + 3.94,
    c = -1/r,
    cnpl = -1/rnpl, 
    cguest = ((cnpl  - (1-frac_guest.x)*(-0.04075193))/frac_guest.x)
  ) #%>%
  select(lipid_guest, cguest)


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

library(broom)  # for tidy model summaries

# Filter data
filtered <- spacings %>%
  filter(temp == 35)

# Fit linear models grouped by lipid_guest and ts
lm_results <- filtered %>%
  filter(lipid_guest %in% c("DOPE", "DOPI (10%)", "Soy PI (10%)", "DOPC (20%)")) %>%
  group_by(lipid_guest, ts) %>%
  do(tidy(lm(cguest ~ press, data = .))) %>%
  filter(term == "(Intercept)")

ggplot(
  lm_results, aes(x = lipid_guest, y = estimate)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  geom_errorbar(aes(ymin = estimate - std.error,
                    ymax = estimate + std.error),
                width = 0.2, color = "black") +
  labs(x = "Lipid Guest", y = expression(paste(c[0], " (1/Å)"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


############################################################################
#Soy PI phase
#Supp Fig 2
#######################################################################################
################################################################################
####### Pt1 is phase, pt2 is dspace + c0 ############################
###############################################################
files_profs = list.files(here("PI-profiles"), pattern = "*.dat", full.names = TRUE) #PI pressure data
file_metdata = here("metadata.csv")

# factor order for guest lipids
lipids_ord = c("Soy PI")


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
  filter(ts == 0) %>%
  #filter(temp == 35) %>%
  #filter(press == 0) %>%
  #filter(pdir == "up") %>%
  #filter(!(samp %in% paste("DM", str_pad(seq(11, 15), 4, pad='0'), sep=''))) %>% 
  # crop q-domain for normalization
  filter(between(q, 0.02, 0.35)) %>% 
  # normalize I(q)
  group_by(samp, press, temp) %>% 
  mutate(
    iq = iq - min(iq),
    iq = iq/max(iq),
    iq = iq + miniq
  )

# average replicate shots
saxsmean = saxsnorm %>% 
  #filter(press == 0) %>% 
  #filter(pdir == "up") %>%
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


# Normalize only within each temperature group
saxsnorm_temp <- saxsmeta %>% 
  filter(ts == 0) %>%
  filter(between(q, 0.02, 0.35)) %>%
  group_by(temp) %>%                              # group by temperature ONLY
  mutate(
    iq = iq - min(iq, na.rm = TRUE),              # shift to zero baseline
    iq = iq / max(iq, na.rm = TRUE),              # scale so max per temp = 1
    iq = iq + 0.1                                 # miniq offset if needed
  ) %>%
  ungroup()

# Average replicate shots
saxsmean_temp <- saxsnorm_temp %>% 
  mutate(rep = as.numeric(rep)) %>% 
  group_by(samp, lipid_guest, frac_guest, ts, temp, press, pdir, q) %>% 
  filter(rep > (max(rep) - 3)) %>% 
  summarise(
    iq = mean(iq),
    err = sqrt(sum(err**2)) / sqrt(n()),
    fname = list(fname),
    .groups = "drop"
  ) %>% 
  mutate(
    lipid_guest = factor(lipid_guest, levels = lipids_ord)
  )

# Plot normalized by temperature
saxsmean_temp %>%
  filter(
    press < 101,
    !(press %% 100),
    pdir == "up",
    lipid_guest != "NA"
  ) %>%
  ggplot(aes(x = q, y = iq)) +
  facet_grid(rows = vars(press*10), cols = vars(temp)) +
  scale_y_log10() +
  geom_line(aes(group = samp)) +
  labs(
    x = "q (1/Å)",
    y = "I(q) (normalized within temp)"
  ) +
  xlim(0.0, 0.35) +
  theme_pubr(x.text.angle = 45) +
  theme_classic(base_size = 6) +
  theme(legend.position = "none")


