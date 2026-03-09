# plotting helper for the "cumulative stacked barplot"
# used to visualize curvature contributions

# Some color palettes riff on the ones in here
source(here("00-scripts", "seelipids_helpers.R")) 

# this gets the collected lit values from the Gsheet etc.
# for reasons I don't totally understand, I have to go into the first script
# and run the last block directly to get c0_allschemes
source(here("00-scripts", "c0_meta_analysis.R")) 
source(here("00-scripts", "Tm_meta_analysis.R")) 

# factor order for headgroups based on curvature
class_order_c0 = c0_allschemes %>% 
  # seems most intuitive
  filter(scheme == "linreg") %>% 
  filter(
    (carbsn1 %in% c(18, NA)),
    (dbonsn1 %in% c(1, NA))
  ) %>% 
  arrange(-c0) %>% 
  distinct(class) %>% 
  .$class

# function to calc contributions to PLCI (PhosphoLipid Curvature Index)
# at present, this requires full chain-resolution
# requires cols: class, carbsn1, dbonsn1
# adds cols: plci, ctol. ctol is measurement error on the PLCI contribution
calc_plci = function(
    # phospholipid structure and abundance data
    pldata,
    # sn-1 chain factors to use
    scheme = c("oleoyl", "satsn1", "linreg")
){
  ## remove sn1 chain data for anionic PLs bc they are not predictors in the model
  #pldata_anionna = pldata %>% mutate(
  #  carbsn1 = ifelse(class %in% c("PG", "PS"), NA, carbsn1),
  #  dbonsn1 = ifelse(class %in% c("PG", "PS"), NA, dbonsn1)
  #)
  
  plcidata = pldata %>% 
    # remove chain data for classes that are not part of the linreg
    mutate(
      carbsn1 = ifelse(!(class %in% class_linreg), NA, carbsn1),
      dbonsn1 = ifelse(!(class %in% class_linreg), NA, dbonsn1)
    ) %>% 
    # duplicate the input data for each scheme to be evaluated
    cross_join(
      tibble(scheme = scheme)
    ) %>% 
    # join the headgroup-only schemes
    left_join(
      c0_allschemes %>% 
        filter(scheme %in% c("oleoyl", "satsn1")) %>% 
        select(scheme, class, c0, tol),
      by = c("scheme", "class")
    ) %>% 
    # join sn1-aware scheme
    mutate(unsatsn1 = as.logical(dbonsn1)) %>% 
    left_join(
      c0_allschemes %>% 
        filter(scheme == "linreg") %>% 
        # remove chain data for classes that are not part of the linreg
        # needed so that the classes still match up!
        mutate(
          carbsn1 = ifelse(!(class %in% class_linreg), NA, carbsn1),
          dbonsn1 = ifelse(!(class %in% class_linreg), NA, dbonsn1)
        ) %>% 
        select(scheme, class, carbsn1, dbonsn1, c0, tol) %>% 
        mutate(unsatsn1 = as.logical(dbonsn1)) %>% 
        select(-dbonsn1),
      by = c("scheme", "class", "carbsn1", "unsatsn1"),
      suffix = c('', ".linreg")
    ) %>% 
    mutate(
      c0  = ifelse(is.na(c0 ), c0.linreg, c0 ),
      tol = ifelse(is.na(tol), tol.linreg,  tol)
    ) %>% 
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
  
  return(plcidata)
}

# custom "cumulative stacked barplot" wrapper
# for visualizing curvature contributions
gg_plcurv = function(
    data,
    darkmode = FALSE,
    thres_draw = 0, # can pass a mole fraction threshold below which color blocks are removed.
    label_frac = 0.015, # min fraction at which bar get labeled; set to 1 to disable labels.
    label_size = 1.5,
    baseline = 0.5/.pt, # *weight* of the baseline
    # aesthetics get passed in as naked args
    # requires `cols = [sp|eid]`
    ...
){
  # unpack the ellipsis args as strings
  mapstrs = lapply(rlang::enexprs(...), as.character)
  # set up row faceting if desired
  if("rows" %in% names(mapstrs)){
    rowvar = vars(eval(sym(mapstrs$rows)))
    rownam = sym(mapstrs$rows)
  }else{
    rowvar = NULL
    rownam = NULL
  }
  # ensure ordering
  data = data %>% arrange(class) # and maybe even id?
  this_gg = data %>%
    # apply threshold
    filter(abs(eval(sym(mapstrs$y))) >= thres_draw) %>% 
    # group by the passed x aesthetic
    #group_by(eval(sym(mapstrs$cols))) %>%
    group_by(eval(sym(mapstrs$cols)), eval(rownam)) %>%
    # sum or mean?
    #summarize(plci = sum(plci)) %>% # try leaving this for the user to do upstream
    mutate(neg = (plci<0)*0.3) %>% # 0.3 is the offset width
    arrange(neg, class) %>% 
    mutate(
      # running total
      end = cumsum(plci),
      y = end - plci/2
    ) %>% 
    ggplot(aes(x = neg)) +
    # due to plotting mechanics, columns are mapped as cols, not x
    facet_grid(
      rows = rowvar,
      cols = vars(eval(sym(mapstrs$cols))), 
      switch = 'x'
    ) +
    geom_hline(
      size = baseline,
      yintercept = 0, color = ifelse(darkmode, "white", "black")
    ) +
    geom_tile(
      aes(
        y = y,
        height = plci,
        fill = factor(class, levels = c("LPC", "LPE", "PS", "PC", "PG", "PI", "PE", "O-PC", "P-PC", "P-PE", "O-PE")),
        #color = class
      ),
      width = 0.25,
      size = 0.05,
      color = ifelse(darkmode, "black", "white")
    ) +
    scale_fill_manual(values = chroma_cl) +
    #scale_color_manual(values = c(TRUE = chroma_ol_dark, FALSE = chroma_ol_lite)[darkmode]) +
    ## ind'lwise error bars
    #geom_errorbar(
    #  data = data %>% 
    #    group_by(sp, sp_eid, class) %>%
    #    summarize(frac_molar = sum(frac_molar)) %>%
    #    full_join(icurv, by = c("class")) %>% 
    #    mutate(plci = frac_molar * c0) %>% 
    #    drop_na() %>% 
    #    group_by(sp, sp_eid, scenario) %>% 
    #    summarize(plci = sum(plci)) %>% 
    #    # average by species
    #    group_by(scenario, sp) %>% 
    #    summarize(
    #      plci_mean = mean(plci),
    #      plci_serr = sd(plci)/sqrt(n()),
    #      neg = TRUE * 0.3
    #    ),
    #  aes(
    #    ymin = plci_mean - plci_serr,
    #    ymax = plci_mean + plci_serr
    #  ),
    #  width = 0.05
    #) +
    theme_pubr() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    guides(color = "none") +
    labs(y = "Mean c0 (1/Å)")
  this_gg
}

# function to calc contributions to PLFI (PhosphoLipid Fluidity Index)
# at present, this requires full chain-resolution,
# though it can be inferred before piping in
# requires cols: class, carbsn1, dbonsn1, carbsn2, dbonsn2
calc_plfi = function(
    # phospholipid structure and abundance data
  pldata,
  # sn-1 chain factors to use
  schemes = "logcanonical"
){
  pldata %>% 
    # take care of the underspecification problem for double bonds
    # by capping at the max number in the training set
    rowwise() %>% 
    mutate(
      dbonsn1_orig = dbonsn1,
      dbonsn2_orig = dbonsn2,
      dbonsn1 = min(dbonsn1, max(data_melt$dbonsn1)),
      dbonsn2 = min(dbonsn2, max(data_melt$dbonsn2)),
      across(contains("dbon"), as.factor),
    ) %>% 
    ungroup() %>% 
    cross_join(
      mods_plfi %>% 
        filter(scheme %in% schemes) %>% 
        select(mod)
    ) %>% 
    group_by(scheme) %>% 
    mutate(
      # keep all the unk Tm lipids, put NAs
      #tm = safely(
      #  predict, 
      #  otherwise = NA
      #)(
      #  mod, 
      #  newdata = cur_data()
      #) %>% print() %>% .$result,
      # the model wants class to be character, for versatility
      tm = predict(
        mod[[1]], 
        # this is how we deal with classes not in the model
        newdata = cur_data() %>% 
          mutate(
            class = as.character(class),
            class = ifelse(class %in% data_melt$class, class, NA)
          )# %>% 
          #print() %>% View()
      ),
      plfi = tm * frac_molar
    )
}

# list of headgroups in fluidity order from least to most
ord_tm = crossing(
  class = names(chroma_cl),
  carbsn1 = 12,
  carbsn2 = 12,
  dbonsn1 = 0,
  dbonsn2 = 0,
  frac_molar = 0
) %>% 
  calc_plfi() %>% 
  drop_na() %>% 
  arrange(-tm) %>% 
  .$class

# a visualization of PLFI more in the style of gg_acylch()
# the problem with it is that trends in the mean are hard to see
# because the x axis is so broad
gg_plvisc = function(data, binwidth = 5, ...){
  # extract the full range for axis setup
  tm = data$tm
  # precalculate the Tm bin breaks
  tm_bks = seq(
    binwidth*floor(min(tm, na.rm = TRUE)/binwidth), 
    binwidth*ceiling(max(tm, na.rm = TRUE)/binwidth), 
    binwidth
  )
  # and the bin midpoints
  tm_mid = (2*tm_bks + binwidth)/2
  # remove last element
  tm_mid = head(tm_mid, -1)
  
  # bin the Tms as specified
  data_histo = data %>% 
    mutate(
      tm_bin = cut(
        tm, 
        breaks = tm_bks,
        labels = tm_mid
      ) %>% 
        as.character() %>% 
        as.numeric()
    ) %>% 
    # and sum up each class in each bin
    group_by(tm_bin, class, .add = TRUE) %>% 
    summarise(frac_molar = sum(frac_molar)) %>% 
    # then norm to relative frequency
    ungroup(class) %>% print() %>% 
    mutate(frac_molar_tot = sum(frac_molar)) %>% 
    ungroup(tm_bin) %>% print() %>% 
    mutate(frac_molar = frac_molar/max(frac_molar_tot))
  
  data_histo %>% 
    gg_acylch(
      cols = NA,
      breakfun = waiver(),
      ...
    )
}

# custom "cumulative stacked barplot" wrapper
# for visualizing PLFI contributions
gg_plfi = function(
    data,
    darkmode = FALSE,
    thres_draw = 0, # can pass a mole fraction threshold below which color blocks are removed.
    label_frac = 0.015, # min fraction at which bar get labeled; set to 1 to disable labels.
    label_size = 1.5,
    baseline = 0.5/.pt, # *weight* of the baseline
    posfirst = TRUE, # show positive contributions coming up from baseline
    # aesthetics get passed in as naked args
    # requires `cols = [sp|eid]`
    ...
){
  # unpack the ellipsis args as strings
  mapstrs = lapply(rlang::enexprs(...), as.character)
  # set up row faceting if desired
  if("rows" %in% names(mapstrs)){
    rowvar = vars(eval(sym(mapstrs$rows)))
    rownam = sym(mapstrs$rows)
  }else{
    rowvar = NULL
    rownam = NULL
  }
  # ensure ordering
  data = data %>% arrange(class) # and maybe even id?
  this_gg = data %>%
    # apply threshold
    filter(abs(eval(sym(mapstrs$y))) >= thres_draw) %>% 
    # group by the passed x aesthetic
    #group_by(eval(sym(mapstrs$cols))) %>%
    group_by(eval(sym(mapstrs$cols)), eval(rownam)) %>%
    # sum or mean?
    mutate(neg = (posfirst-0.5)*2*(plfi<0)*0.3) %>% # 0.3 is the offset width
    arrange(neg, class) %>% 
    mutate(
      # running total
      end = cumsum(plfi),
      y = end - plfi/2
    ) %>% 
    ggplot(aes(x = neg)) +
    # due to plotting mechanics, columns are mapped as cols, not x
    facet_grid(
      rows = rowvar,
      cols = vars(eval(sym(mapstrs$cols))), 
      switch = 'x'
    ) +
    geom_hline(
      size = baseline,
      yintercept = 0, color = ifelse(darkmode, "white", "black")
    ) +
    geom_tile(
      aes(
        y = y,
        height = plfi,
        fill = class,
        #color = class
      ),
      width = 0.25,
      size = 0.05,
      color = ifelse(darkmode, "black", "white")
    ) +
    scale_fill_manual(values = chroma_cl) +
  theme_pubr() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    guides(color = "none") +
    labs(y = "PLFI (°C)")
  this_gg
}
