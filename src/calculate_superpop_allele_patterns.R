# calculate_allele_patterns.R

# script to take allele counts per locus
# and calculate the probability of the observed pattern

# Daniel Cotter

## load packages -------
suppressPackageStartupMessages(require(tidyverse))
library(optparse)

## parse arguments ------
option_list <- list(
  make_option(c("--chr"),
              type = "character", default = NULL,
              help = "chromosome to analyze"),
  make_option(c("--drop-singletons"),
              type= "logical", action = "store_true", default = FALSE,
              help = "flag specifying whether to drop singletons"),
  make_option(c("--threshold"),
              type = "numeric", default = 0.05,
              help = "threshold to define a rare vs common allele"),
  make_option(c("--sample"),
              type = "numeric", default = 0,
              help = "number of snps to sample. default is whole chromosome")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if ( is.null(opt$chr) ) {
  print_help(opt_parser)
  stop("Must specify a chromosome (e.g. --chr 22).", call.=FALSE)
}

pop_label <- "superpops"
CHR <- opt$chr
DROP_SINGLETONS <- opt$`drop-singletons`
z <- opt$threshold
g_list <- seq(10, 300, by = 10)
sample_size = opt$sample

## create output directory if non existent -----
if (!file.exists(file.path("data", "patterns"))) {
  dir.create(file.path("data", "patterns"), showWarnings = FALSE)
}

## read in superpop data -----
df <- read.table(
  file.path("data", "allele_counts", 
  paste("chr", CHR, "_counts_", pop_label, ".txt", sep="")),
  header = TRUE)

# subset if sample size > 0
if (sample_size > 0) {
  set.seed(42)
  df_sub <- df %>% sample_n(sample_size, replace=FALSE)
  df <- df_sub %>% tibble()
  sample_label <- paste(format(sample_size, scientific=FALSE), '-snps', sep="")
} else {
  sample_label <- "all-snps"
}


df_long <- df %>% 
  gather(pop, counts, -c(chr:major)) %>% 
  select(-c(minor:major)) %>% 
  separate(counts, into = c("minor", "major"), sep = "/") %>%
  mutate(minor = as.integer(minor), major = as.integer(major))

rm(list = c("df"))

# remove alleles for which there is no globally "minor" allele (i.e both = 50%)
df_long <- df_long %>% filter(!is.na(minor))

# drop non-biallelic sites
drop <- df_long %>% 
  group_by(pos) %>% 
  summarise(tot_minor = sum(minor)) %>% 
  filter(tot_minor == 0) %>% 
  pull(pos)
df_long <- df_long %>% filter(!(pos %in% drop))

# drop singletons if necessary
if (DROP_SINGLETONS) {
  drop1 <- df_long %>% 
    group_by(pos) %>% 
    summarise(tot_minor = sum(minor)) %>% 
    filter(tot_minor==1) %>% 
    pull(pos)
  df_long <- df_long %>% filter(!(pos %in% drop1))
}

# ## grab the actual patterns for each locus
actual_patterns <- df_long %>%
  mutate(freq = minor / (minor + major)) %>%
  mutate(code = ifelse(freq == 0, "U", ifelse(freq <= 0.05, "R", "C"))) %>%
  select(-c(tot_alleles, minor, major, freq)) %>%
  spread(pop, code) %>%
  mutate(pattern=paste(AFR, EUR, SAS, EAS, AMR, sep = "")) %>%
  select(pos, pattern)

actual_pattern_vector <- actual_patterns %>%
  group_by(pattern) %>% 
  summarise(n=n()) %>% 
  mutate(freq=n/sum(n)) %>% 
  arrange(-freq)

## write the actual pattern probabilities to file
if (DROP_SINGLETONS) {
  write.table(actual_pattern_vector, file = paste('data/patterns/',
                                     CHR, '_actualPatternn_', sample_label, 
                                     '_noSingletons.txt', sep=""),
              sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
} else {
  write.table(actual_pattern_vector, file = paste('data/patterns/',
                                     CHR, '_actualPattern_', sample_label, 
                                     '_wSingletons.txt', sep=""),
              sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
}

## define functions to calculate the probability of (U)nobserved or (R)are -----
# (C)ommon is defined by 1 - U - R

u_matrix <- function(minor, major, g) {
  u_prob <- function(minor, major, g) {
    g <- rep(g, length(minor))
    U <- choose(major, g) / choose(minor + major, g)
    return(U)
  }
  mat <- data.frame(minor = minor, major = major) %>%
  dplyr::distinct() %>%
  mutate(U = round(u_prob(minor, major, g), 6)) %>%
  tidyr::pivot_wider(names_from = major, values_from = U) %>%
  tibble::column_to_rownames("minor")
  return(as.matrix(mat))
}

r_matrix <- function(minor, major, g, z) {
  r_prob <- function(minor, major, g, z) {
    g <- rep(g, length(minor))
    cutoff <- floor(z * (minor + major))
    R <- mapply(function(n1, n2, g, k) {
      i <- 1:k
      sum(choose(n1, i) * choose(n2, g - i))
      },
      minor, major, g, cutoff,
      SIMPLIFY = TRUE) / (choose(minor + major, g))
    return(R)
  }
  mat <- data.frame(minor = minor, major = major) %>%
  dplyr::distinct() %>%
  mutate(R = round(r_prob(minor, major, g, z), 6)) %>%
  tidyr::pivot_wider(names_from = major, values_from = R) %>%
  tibble::column_to_rownames("minor")
  return(as.matrix(mat))
}

## define function to calculate probabilites   ----
merge_codes <- function(..., log = FALSE, patterns) {
  x <- data.frame(...) %>% expand.grid()
  if (log) {
    products <- log(x[[1]])
    for (i in 2:length(x)) {
      products <- products + log(x[[i]])
    }
  } else {
    products <- x[[1]]
    for (i in 2:length(x)) {
      products <- products * x[[i]]
    }
  }
  return(data.frame(pattern = patterns, prob = products)) 
  #%>% filter(prob > 0))
}

## calculate the three probabilites ----------
df_all <- NULL
df_match <- NULL

for (i in seq_along(g_list)) {
  g <- g_list[i]
  print(g)
  # calculate lookup matrices for R and U
  U_mat <- u_matrix(df_long$minor, df_long$major, g)
  R_mat <- r_matrix(df_long$minor, df_long$major, g, z)

  # use the lookup matrices to calculat U, R, and C 
  # for each position / population pair
  df_probs <- df_long %>%
  mutate(
    U = U_mat[cbind(as.character(minor), as.character(major))],
    R = R_mat[cbind(as.character(minor), as.character(major))]) %>%
  mutate(C=round(1-R-U, 6)) %>%
  mutate(C = sapply(C, function(x) { ifelse(x < 0, 0, x) })) %>%
  select(-minor, -major) %>%
  gather(cat, prob, R, U, C) %>%
  spread(pop, prob)
  # filter for loci where N_j for any population is less than g and remove them
  # (on chr22 only 115 loci have <300 for any given pop)
  df_probs <- df_probs %>% filter(across(everything(), ~!is.na(.)))

  # get column of patterns
  pattern_vec <- data.frame(a = c("C", "R", "U"),
                            b = c("C", "R", "U"),
                            c = c("C", "R", "U"),
                            d = c("C", "R", "U"),
                            e = c("C", "R", "U")) %>%
                            expand.grid() %>%
    unite(pattern, a:e, sep="") %>%
    pull(pattern)

  df_probs <- df_probs %>%
  group_by(chr, pos, tot_alleles) %>%
  arrange(pos, cat) %>%
    do(merge_codes(.$AFR, .$EUR, .$SAS, .$EAS, .$AMR,
     log = FALSE, patterns = pattern_vec))

  if (DROP_SINGLETONS) {
    write.table(df_probs, file = paste('data/patterns/',
                                       CHR, '_', g, '_pattern_byPosition_', sample_label, 
                                       '_noSingletons.txt', sep=""),
                sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
  } else {
    write.table(df_probs, file = paste('data/patterns/',
                                       CHR, '_', g, '_pattern_byPosition_', sample_label, 
                                       '_wSingletons.txt', sep=""),
                sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
  }

  df_patterns <- df_probs %>%
    spread(pattern, prob) %>%
    group_by(chr) %>%
    summarise(across(everything(), ~mean(.x))) %>%
    gather(pattern, prob, -c(chr, pos, tot_alleles)) %>%
    select(-c(chr, pos, tot_alleles)) %>%
    rename(!!paste(g) := prob)

  if (is.null(df_all)) {
    df_all <- df_patterns
  } else {
    df_all <- inner_join(df_all, df_patterns, by="pattern")
  }

  ## check if the most likely prob (sans UUUUU) matches the actual probability
  df_rate <- df_probs %>%
    group_by(pos) %>%
    filter(pattern!='UUUUU') %>%
    filter(prob==max(prob)) %>%
    rename(pattern2=pattern) %>%
    inner_join(actual_patterns, by="pos") %>%
    mutate(match=(pattern==pattern2)) %>%
    summarise(match=sum(match)) %>%
    ungroup() %>%
    summarise(match_rate=mean(match)) %>%
    mutate(g=paste(g))

  if (is.null(df_match)) {
    df_match <- df_rate
  } else {
    df_match <- rbind(df_match, df_rate)
  }
}

if (DROP_SINGLETONS) {
  write.table(df_all, file = paste('data/patterns/',
                                   CHR, '_patterns_', sample_label, '_noSingletons.txt', sep=""),
              sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
} else {
  write.table(df_all, file = paste('data/patterns/',
                                     CHR, '_patterns_', sample_label, '_wSingletons.txt', sep=""),
              sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
}

if (DROP_SINGLETONS) {
  write.table(df_match, file = paste('data/patterns/',
                                   CHR, '_pattern-match-proportions_', sample_label, '_noSingletons.txt', sep=""),
              sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
} else {
  write.table(df_match, file = paste('data/patterns/',
                                   CHR, '_pattern-match-proportions_', sample_label, '_wSingletons.txt', sep=""),
              sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
}
