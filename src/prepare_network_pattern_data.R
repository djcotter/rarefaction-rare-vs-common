## load packages ------
library(tidyverse)
library(optparse)

## parse arguments ------
option_list <- list(
  make_option(
    c("--input"),
    type = "character",
    default = NULL,
    help = "path to input file"
  ),
  make_option(
    c("--output"),
    type = "character",
    default = NULL,
    help = "path to output file"
  ),
  make_option(
    c("--relative"),
    type = "logical",
    action = "store_true",
    default = FALSE,
    help = "flag specifying whether to normalize the data"
  ),
  make_option(
    c("--g_filter"),
    type = "numeric",
    default = NULL,
    help = "g column to filter out and group"
  )
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$wSingletons)) {
  print_help(opt_parser)
  stop("Must provide a file with patterns including singletons", call. =
         FALSE)
}

if (is.null(opt$noSingletons)) {
  print_help(opt_parser)
  stop("Must provide a file with patterns without singletons", call. = FALSE)
}

## read in data -------
df <-  read.table(opt$input, header = T, check.names = FALSE)

if (is.null(opt$g_filter)) {
  g_filter <- max(as.numeric(names(df)[2:length(names(df))]))
} else {
  g_filter <- opt$g_filter
}

## define function to combine patterins into ordered triples
combine_pattern <- function(pattern) {
  new_pattern <- paste(str_count(pattern, 'U'),
                       str_count(pattern, 'R'),
                       str_count(pattern, 'C'),
                       sep = ',')
  return(new_pattern)
}

## modify the data and calculate pattern proportions
df_mod <-
  df %>% select(pattern,!!paste(g_filter)) %>%
  rename(node = pattern, prop = !!paste(g_filter)) %>%
  mutate(node = combine_pattern(node)) %>%
  group_by(node) %>%
  summarise(prop = sum(prop))

if (opt$relative) {
  df_mod_relative <-
    df_mod %>% filter(node != '5,0,0') %>%
    mutate(prop = prop / sum(prop))
  df_mod_relative <-
    rbind(df_mod_relative, data.frame(node = '5,0,0', prop = 0))
  write.table(
    df_mod_relative,
    file = opt$output,
    col.names = T,
    row.names = F,
    quote = F,
    sep = '\t'
  )
} else {
  write.table(
    df_mod,
    file = opt$output,
    col.names = T,
    row.names = F,
    quote = F,
    sep = '\t'
  )
}
