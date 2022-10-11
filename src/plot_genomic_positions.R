# calculate_allele_patterns.R
# script to take allele counts per locus and calculate the probability of the observed pattern
# Daniel Cotter

## load packages -------
suppressPackageStartupMessages(require(tidyverse))
library(ggpubr)
library(RColorBrewer)
#suppressPackageStartupMessages(require(optparse))

## temp stuff
setwd('~/Projects/rarefaction-project/')
pop_label = "superpops"
CHR = 1
g = 300
z = 0.05
DROP_SINGLETONS = TRUE

## read in superpop data -----
df <- read.table(file.path('data', 'allele_counts', paste('chr', CHR, '_counts_', pop_label, '.txt', sep="")), header = TRUE)

# temporarily subset the sample size down
# set.seed(1)
# df_sub <- df %>% sample_n(100000, replace=FALSE)
# df <- df_sub %>% tibble()

df_long <- df %>% 
  gather(pop, counts, -c(chr:major)) %>% 
  select(-c(minor:major)) %>% 
  separate(counts, into=c('minor', 'major'), sep="/") %>%
  mutate(minor=as.integer(minor), major=as.integer(major))

rm(list=c('df'))

# remove alleles for which there is no globally "minor" allele (i.e both = 50%)
df_long <- df_long %>% filter(!is.na(minor))

# drop non-biallelic sites
drop <- df_long %>% 
  group_by(pos) %>% 
  summarise(tot_minor = sum(minor)) %>% 
  filter(tot_minor==0) %>% 
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
# actual_patterns <- df_long %>% 
#   mutate(freq=minor/(minor+major)) %>% 
#   mutate(code=ifelse(freq==0, 'U', ifelse(freq<0.05, 'R', 'C'))) %>% 
#   select(-c(tot_alleles, minor, major, freq)) %>% 
#   spread(pop, code) %>% 
#   mutate(pattern=paste(AFR, EUR, SAS, EAS, AMR, sep="")) %>% 
#   select(pos, pattern)

## define functions to calculate the probability of (U)nobserved or (R)are -----
# (C)ommon is defined by 1 - U - R

U_matrix <- function(minor, major, g) {
  U_prob <- function(minor, major, g) {
    g = rep(g, length(minor))
    U <- choose(major, g) / choose(minor+major, g)
    return(U)
  }
  mat <- data.frame(minor=minor,major=major) %>% distinct()
  mat <- mat %>% mutate(U=round(U_prob(minor, major, g),6))
  mat <- mat %>% pivot_wider(names_from = major, values_from = U)
  mat <- mat %>% column_to_rownames('minor')
  return(as.matrix(mat))
}

R_matrix <- function(minor, major, g, z) {
  R_prob <- function(minor, major, g, z) {
    g <- rep(g, length(minor))
    cutoff <- ceiling(z*(minor+major)) -1
    R <- mapply(function(N1, N2, g, k) {i=1:k; sum( choose(N1, i) * choose(N2, g-i) )}, minor, major, g, cutoff, SIMPLIFY = TRUE) / (choose(minor+major, g))
    return(R)
  }
  mat <- data.frame(minor=minor,major=major) %>% distinct()
  mat <- mat %>% mutate(R=round(R_prob(minor, major, g, z),6))
  mat <- mat %>% pivot_wider(names_from = major, values_from = R)
  mat <- mat %>% column_to_rownames('minor')
  return(as.matrix(mat))
}

## define function to calculate probabilites   ----
merge_codes <- function(..., log=FALSE, patterns) {
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
  return(data.frame(pattern=patterns,prob=products) %>% filter(prob>0))
}

## calculate the three probabilites ----------
print(g)
# calculate lookup matrices for R and U
U_mat <- U_matrix(df_long$minor, df_long$major, g)
R_mat <- R_matrix(df_long$minor, df_long$major, g, z)

# use the lookup matrices to calculat U, R, and C for each position / population pair
df_probs <-
  df_long %>% mutate(U = U_mat[cbind(as.character(minor), as.character(major))],
                     R = R_mat[cbind(as.character(minor), as.character(major))]) %>%
  mutate(C = round(1 - R - U, 6)) %>%
  mutate(C = sapply(C, function(x) {
    ifelse(x < 0, 0, x)
  })) %>%
  select(-minor,-major) %>%
  gather(cat, prob, R, U, C) %>%
  spread(pop, prob)
# filter for loci where N_j for any population is less than g and remove them
# (on chr22 only 115 loci have <300 for any given pop)
df_probs <- df_probs %>% filter(across(everything(), ~ !is.na(.)))

# get column of patterns
pattern_vec <- data.frame(
  a = c("C", "R", "U"),
  b = c("C", "R", "U"),
  c = c("C", "R", "U"),
  d = c("C", "R", "U"),
  e = c("C", "R", "U")
) %>% expand.grid() %>%
  unite(pattern, a:e, sep = "") %>% pull(pattern)


df_probs <-
  df_probs %>% group_by(chr, pos, tot_alleles) %>% arrange(pos, cat) %>%
  do(merge_codes(
    .$AFR,
    .$EUR,
    .$SAS,
    .$EAS,
    .$AMR,
    log = FALSE,
    patterns = pattern_vec
  ))

 if (DROP_SINGLETONS) {
  write.table(
    df_probs,
    file = paste(
      '~/Projects/rarefaction-project/data/tmp/',
      CHR,
      '_',
      g,
      '_pattern_byPosition_noSingletons.txt',
      sep = ""
    ),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
} else {
  write.table(
    df_probs,
    file = paste(
      '~/Projects/rarefaction-project/data/tmp/',
      CHR,
      '_',
      g,
      '_pattern_byPosition.txt',
      sep = ""
    ),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

combine_pattern <- function(pattern) {
  new_pattern <- paste(
    str_count(pattern, 'C'),
    'C',
    str_count(pattern, 'R'),
    'R',
    str_count(pattern, 'U'),
    'U',
    sep = ''
  )
  return(new_pattern)
}

## g = 300 (w/ singletons) ---------------------
chr_patterns <- read.table(paste(
  '~/Projects/rarefaction-project/data/tmp/',CHR,'_', g, '_pattern_byPosition.txt', sep=""),
  header = T)
# set the size of the windows here
win_size_txt = '100kb'
win_size=100000
windows <- (plyr::round_any(chr_patterns$pos - (min(chr_patterns$pos)-1),win_size,f=ceiling) - win_size/2)/1e6
chr_patterns$windows <- windows

## Combine Patterns - w/ UUUUU
chr_window_pattern_vec <- chr_patterns %>% 
  dplyr::group_by(windows) %>% 
  mutate(win_sum = sum(prob)) %>% 
  ungroup() %>% group_by(windows, pattern) %>% 
  mutate(pattern_sum = sum(prob)) %>% 
  ungroup() %>% mutate(win_prob = pattern_sum/ win_sum) %>% 
  distinct(pattern, windows, .keep_all = T) %>% 
  select(-tot_alleles, -pattern_sum, -prob, -win_sum) %>% 
  mutate(pattern=combine_pattern(pattern))

# color_pallete
color_vec <- colorRampPalette(brewer.pal(12,"Set3"))(length(unique(chr_window_pattern_vec$pattern)))

p1 <- ggplot(chr_window_pattern_vec %>% 
               mutate(pattern=fct_rev(pattern)), 
             aes(x=windows, y=win_prob, fill=pattern)) + 
  geom_col() + scale_fill_manual('Pattern', values=color_vec) +6 
  theme_pubr() + xlab('Position (Mb)') + ylab('Average Probability') +
  guides(fill=guide_legend(ncol=1))

## Combine Patterns no UUUUU
chr_window_pattern_vec_noU <- chr_patterns %>% 
  filter(pattern!='UUUUU') %>%
  dplyr::group_by(windows) %>% 
  mutate(win_sum = sum(prob)) %>% 
  ungroup() %>% group_by(windows, pattern) %>% 
  mutate(pattern_sum = sum(prob)) %>% 
  ungroup() %>% mutate(win_prob = pattern_sum/ win_sum) %>% 
  distinct(pattern, windows, .keep_all = T) %>% 
  select(-tot_alleles, -pattern_sum, -prob, -win_sum) %>%
  mutate(pattern=combine_pattern(pattern))

# color_pallete
color_vec <- color_vec[1:length(color_vec)-1]

p2 <- ggplot(chr_window_pattern_vec_noU %>% 
               mutate(pattern=fct_rev(pattern)), 
             aes(x=windows, y=win_prob, fill=pattern)) + 
  geom_col() + scale_fill_manual('Pattern', values=color_vec) +
  theme_pubr() + xlab('Position (Mb)') + ylab('Relative Probability') + 
  guides(fill=guide_legend(ncol=1))

## Max patterns (ungrouped)
chr_max_patterns <- chr_patterns %>% 
  filter(pattern!='UUUUU') %>% 
  group_by(windows) %>% 
  count(pattern) %>% 
  top_n(1)

p3 <- ggplot(chr_max_patterns,
             aes(x=windows,
                 y=fct_rev(pattern))) + 
  geom_line(aes(group=1)) +
  theme_pubr() + xlab('Position (Mb)') +
  ylab('Most Frequent Pattern')

## Max patterns (grouped)
chr_max_patterns_grouped <- chr_patterns %>% 
  filter(pattern!='UUUUU') %>% 
  group_by(windows) %>% 
  count(pattern) %>% 
  top_n(1) %>%
  mutate(pattern=combine_pattern(pattern))

p4 <- ggplot(chr_max_patterns_grouped,
             aes(x=windows,
                 y=pattern)) + 
  geom_line(aes(group=1)) +
  theme_pubr() + xlab('Position (Mb)') +
  ylab('Most Frequent Pattern')

p <- ggarrange(p1,p2,p4, ncol=1, common.legend = T, legend='right', align = 'v')
ggsave(p,filename = paste('~/../Downloads/chr', CHR, '_patterns_', win_size_txt, '_g', g, '.pdf', sep=""), width = 7, height=7, units='in')


## g = 300 (no Singletons) ----------------------------------------------------
chr_patterns <- read.table(paste(
  '~/Projects/rarefaction-project/data/tmp/',CHR,'_', g, '_pattern_byPosition_noSingletons.txt', sep=""),
  header = T)
# set the size of the windows here
win_size_txt = '100kb'
win_size=100000
windows <- (plyr::round_any(chr_patterns$pos - (min(chr_patterns$pos)-1),win_size,f=ceiling) - win_size/2)/1e6
chr_patterns$windows <- windows

## Combine Patterns - w/ UUUUU
chr_window_pattern_vec <- chr_patterns %>% 
  dplyr::group_by(windows) %>% 
  mutate(win_sum = sum(prob)) %>% 
  ungroup() %>% group_by(windows, pattern) %>% 
  mutate(pattern_sum = sum(prob)) %>% 
  ungroup() %>% mutate(win_prob = pattern_sum/ win_sum) %>% 
  distinct(pattern, windows, .keep_all = T) %>% 
  select(-tot_alleles, -pattern_sum, -prob, -win_sum) %>% 
  mutate(pattern=combine_pattern(pattern))

# color_pallete
color_vec <- colorRampPalette(brewer.pal(12,"Set3"))(length(unique(chr_window_pattern_vec$pattern)))

p1 <- ggplot(chr_window_pattern_vec %>% 
               mutate(pattern=fct_rev(pattern)), 
             aes(x=windows, y=win_prob, fill=pattern)) + 
  geom_col() + scale_fill_manual('Pattern', values=color_vec) +
  theme_pubr() + xlab('Position (Mb)') + ylab('Average Probability') +
  guides(fill=guide_legend(ncol=1))

## Combine Patterns no UUUUU
chr_window_pattern_vec_noU <- chr_patterns %>% 
  filter(pattern!='UUUUU') %>%
  dplyr::group_by(windows) %>% 
  mutate(win_sum = sum(prob)) %>% 
  ungroup() %>% group_by(windows, pattern) %>% 
  mutate(pattern_sum = sum(prob)) %>% 
  ungroup() %>% mutate(win_prob = pattern_sum/ win_sum) %>% 
  distinct(pattern, windows, .keep_all = T) %>% 
  select(-tot_alleles, -pattern_sum, -prob, -win_sum) %>%
  mutate(pattern=combine_pattern(pattern))

# color_pallete
color_vec <- color_vec[1:length(color_vec)-1]

p2 <- ggplot(chr_window_pattern_vec_noU %>% 
               mutate(pattern=fct_rev(pattern)), 
             aes(x=windows, y=win_prob, fill=pattern)) + 
  geom_col() + scale_fill_manual('Pattern', values=color_vec) +
  theme_pubr() + xlab('Position (Mb)') + ylab('Relative Probability') + 
  guides(fill=guide_legend(ncol=1))

## Max patterns (ungrouped)
chr_max_patterns <- chr_patterns %>% 
  filter(pattern!='UUUUU') %>% 
  group_by(windows) %>% 
  count(pattern) %>% 
  top_n(1)

p3 <- ggplot(chr_max_patterns,
             aes(x=windows,
                 y=fct_rev(pattern))) + 
  geom_line(aes(group=1)) +
  theme_pubr() + xlab('Position (Mb)') +
  ylab('Most Frequent Pattern')

## Max patterns (grouped)
chr_max_patterns_grouped <- chr_patterns %>% 
  filter(pattern!='UUUUU') %>% 
  group_by(windows) %>% 
  count(pattern) %>% 
  top_n(1) %>%
  mutate(pattern=combine_pattern(pattern))

p4 <- ggplot(chr_max_patterns_grouped,
             aes(x=windows,
                 y=pattern)) + 
  geom_line(aes(group=1)) + 
  theme_pubr() + xlab('Position (Mb)') +
  ylab('Most Frequent Pattern')

p <- ggarrange(p1,p2,p4, ncol=1, common.legend = T, legend='right', align = 'v')
ggsave(p,filename = paste('~/../Downloads/chr', CHR, '_patterns_', win_size_txt, '_g', g, '_noSingeltons.pdf', sep=""), width = 7, height=7, units='in')



## -------------------------
CHR = 6
g = 300
z = 0.05
DROP_SINGLETONS = TRUE
chr_patterns <- read.table(paste(
  '~/Projects/rarefaction-project/data/tmp/',CHR,'_', g, '_pattern_byPosition_noSingletons.txt', sep=""),
  header = T)
# set the size of the windows here
win_size_txt = '250kb'
win_size=250000
windows <- (plyr::round_any(chr_patterns$pos - (min(chr_patterns$pos)-1),win_size,f=ceiling) - win_size/2)/1e6
chr_patterns$windows <- windows

chr_window_pattern_ranks <- chr_patterns %>% 
  filter(pattern!='UUUUU') %>%
  mutate(pattern=combine_pattern(pattern)) %>%
  dplyr::group_by(windows) %>% 
  mutate(win_sum = sum(prob)) %>% 
  ungroup() %>% group_by(windows, pattern) %>% 
  mutate(pattern_sum = sum(prob)) %>% 
  ungroup() %>% mutate(win_prob = pattern_sum/ win_sum) %>% 
  distinct(pattern, windows, .keep_all = T) %>% 
  select(-tot_alleles, -pattern_sum, -prob, -win_sum) %>% 
  group_by(windows) %>% 
  mutate(pattern_rank = dense_rank(desc(win_prob))) %>% 
  arrange(windows, pattern_rank)


RANK_COLOR_CUTOFF = 3
patterns_to_color = chr_window_pattern_ranks %>% filter(pattern_rank <= RANK_COLOR_CUTOFF) %>% pull(pattern) %>% unique()
all_patterns = chr_window_pattern_ranks %>% pull(pattern) %>% unique()
pattern_vec = c(patterns_to_color, all_patterns[!(all_patterns %in% patterns_to_color)])
color_vec = c(colorRampPalette(brewer.pal(length(patterns_to_color),"Set1"))(length(patterns_to_color)), 
              rep("#D3D3D3", 21-length(patterns_to_color)))
alpha_vec = c(rep(1, length(patterns_to_color)),
              rep(0.075, 21-length(patterns_to_color)))


p_rank <- ggplot(chr_window_pattern_ranks %>%
                   mutate(pattern_rank=fct_rev(as.factor(as.numeric(pattern_rank)))), 
                 aes(x=windows, y=pattern_rank)) + 
  geom_line(aes(group=pattern, color=pattern, alpha=pattern)) +
  scale_color_manual("Pattern", breaks=pattern_vec[1:length(patterns_to_color)], values=color_vec) +
  scale_alpha_manual(,breaks=pattern_vec, values=alpha_vec, guide='none') + 
  xlab('Position (Mb)') +
  ylab('Pattern Rank') +
  theme_pubr(legend = "right")
p_rank 

ggsave(p_rank,filename = paste('~/../Downloads/chr', CHR, '_patterns_', win_size_txt, '_g', g, '_noSingeltons_byRank.pdf', sep=""), 
       width = 140, height=90, units='mm')

