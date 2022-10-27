# plot_patterns_vs_sample_size.R

# script to take pattern probabilities and plot patterns as a function
# of sample size g

# Daniel Cotter

## load packages -------
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(optparse))
library(ggpubr)
library(RColorBrewer)

## load packages -------
suppressPackageStartupMessages(require(tidyverse))
library(optparse)

## parse arguments ------
option_list <- list(
  make_option(c("--wSingletons"),
              type= "character", default = NULL,
              help = "path to patterns calculated INCLUDING singletons"),
  make_option(c("--noSingletons"),
              type= "character", default = NULL,
              help = "path to patterns calculated WITHOUT singletons"),
  make_option(c("--actual-wSingletons"),
              type= "character", default = NULL,
              help = "path to actual empirical patterns WITH singletons"),
  make_option(c("--actual-noSingletons"),
              type= "character", default = NULL,
              help = "path to actual empirical patterns WITHOUT singletons"),
  make_option(c("--prob-threshold"),
              type = "numeric", default = 0.01,
              help = "pattern probability threshold to plot")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if ( is.null(opt$wSingletons) ) {
  print_help(opt_parser)
  stop("Must provide a file with patterns including singletons", call.=FALSE)
}

if ( is.null(opt$noSingletons) ) {
  print_help(opt_parser)
  stop("Must provide a file with patterns without singletons", call.=FALSE)
}

if ( is.null(opt$actual) ) {
  print_help(opt_parser)
  stop("Must provide a file with actual emprical pattern freqs", call.=FALSE)
}

# temp file paths ------ REMOVE 
#opt$wSingletons <- '~/Projects/rarefaction-project/data/patterns/22_patterns_all-snps_wSingletons.txt'
#opt$noSingletons <- '~/Projects/rarefaction-project/data/patterns/22_patterns_all-snps_noSingletons.txt'
opt$`actual-wSingletons` <- '~/Projects/rarefaction-project/data/patterns/22_actualPattern_all-snps_wSingletons.txt'
opt$`actual-noSingletons` <- '~/Projects/rarefaction-project/data/patterns/22_actualPattern_all-snps_noSingletons.txt'

## create output directory if non existent ------------------
if (!file.exists(file.path("data", "figures"))) {
  dir.create(file.path("data", "figures"), showWarnings = FALSE)
}

## read in and format data files ----------------------------
read_and_reformat <- function(filepath, normalize=FALSE) {
  # take a filepath and reorganize the input for plotting
  recolor_patterns <- function(pattern, keep) {
    # change the pattern names to only those in the keep vector
    recolor <- ifelse(pattern %in% keep, pattern, 'Other')
    return(data.frame(recolor=recolor))
  }
  # read in the data frame from file path
  df <- read.table(file=filepath, header = T, check.names = F)
  if (normalize) { # if normalize remove UUUUU and recalculate the other probs
    singletons <- grepl('wSingletons', filepath)
    if (singletons) {
      actual_df <- read.table(file=opt$`actual-wSingletons`, 
                              header = T, check.names = F)
    } else {
      actual_df <- read.table(file=opt$`actual-noSingletons`, 
                              header = T, check.names = F)
    }
    relative_df <- df %>%
      gather(g, prob, -pattern) %>%
      mutate(g=as.numeric(g)) %>%
      spread(g, prob) %>%
      filter(pattern!='UUUUU') %>%
      mutate(across(where(is.numeric), ~./sum(.))) %>%
      gather(g, prob, -pattern)
    keep <- relative_df %>%
      filter(g==300) %>%
      filter(prob >= 0.01) %>%
      pull(pattern)
    relative_df <- relative_df %>%
      do(cbind(., recolor_patterns(.$pattern, keep)))
    keep2 <- relative_df %>% 
      filter(recolor!="Other") %>% 
      pull(recolor) %>% 
      unique()
    max_g <- as.character(relative_df %>% mutate(g=as.numeric(g)) %>% pull(g) %>% unique() %>% max() + 20)
    actual_df <- actual_df %>% select(-n) %>% rename(prob=freq) %>%
      do(cbind(., recolor_patterns(.$pattern, keep2))) %>% 
      mutate(g=max_g)
    df_plot <- rbind(relative_df, actual_df) %>% mutate(g=as.numeric(g))
    return(df_plot)
  } else {
    keep <- df %>% gather(g, prob, -pattern) %>%
      filter(g==300) %>% filter(prob >= 0.01) %>%
      pull(pattern)
    df_plot <- df %>%
      gather(g, prob, -pattern) %>%
      mutate(g=as.numeric(g)) %>%
      do(cbind(., recolor_patterns(.$pattern, keep)))
    return(df_plot)
  }
}

df_wSingletons <- read_and_reformat(opt$wSingletons)
df_wSingletons_relative <- read_and_reformat(opt$wSingletons, normalize = TRUE)
df_noSingletons <- read_and_reformat(opt$noSingletons)
df_noSingletons_relative <- read_and_reformat(opt$noSingletons, normalize = TRUE)

## determine pattern colors ------------------------
all_colored_patterns <- c(df_noSingletons %>% pull(recolor), 
                          df_wSingletons %>% pull(recolor),
                          df_wSingletons_relative %>% pull(recolor),
                          df_noSingletons_relative %>% pull(recolor)) %>% unique()

