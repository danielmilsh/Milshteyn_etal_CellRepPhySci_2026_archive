library(tidyverse)
library(pbapply)

# standard multiread function and a simple 
# positional parser for the filename
read_saxsdats = function(datfiles){
  datfiles %>% 
    read_delim(
      delim = "  ",
      comment = '#',
      col_names = c("q", "iq", "err"),
      col_types = list(col_double(), col_double(), col_double()),
      id = "fname"
    ) %>% 
    # JWL0308_80C_250MPa_up_2_data_000001_00001
    # parse filename
    group_by(fname) %>% 
    mutate(fname = basename(fname)) %>% 
    separate(
      fname,
      into = c(
        "file",
        "samp",
        "temp",
        "press",
        "pdir",
        "rep",
        "typ",
        "ser", # series
        "frm"  # frame
      ),
      remove = FALSE
    ) %>% 
    # get rid of units
    mutate(across(
      all_of(c("temp", "press")),
      function(st){st %>% str_extract("(\\d)+") %>% as.numeric()}
    ))
}
