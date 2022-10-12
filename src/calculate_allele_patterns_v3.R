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
CHR = 22
g_list = seq(10,300,by = 10)
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
df_all <- NULL
df_match <- NULL

for (i in 1:length(g_list)) {
  g = g_list[i]
  print(g)
  # calculate lookup matrices for R and U
  U_mat <- U_matrix(df_long$minor, df_long$major, g)
  R_mat <- R_matrix(df_long$minor, df_long$major, g, z)

  # use the lookup matrices to calculat U, R, and C for each position / population pair
  df_probs <- df_long %>% mutate(U = U_mat[cbind(as.character(minor), as.character(major))], 
                                 R = R_mat[cbind(as.character(minor), as.character(major))]) %>% 
    mutate(C=round(1-R-U, 6)) %>% 
    mutate(C=sapply(C, function(x){ifelse(x<0,0,x)})) %>%
    select(-minor, -major) %>% 
    gather(cat, prob, R, U, C) %>% 
    spread(pop, prob) 
  # filter for loci where N_j for any population is less than g and remove them
  # (on chr22 only 115 loci have <300 for any given pop)
  df_probs <- df_probs %>% filter(across(everything(), ~!is.na(.)))
  
  # get column of patterns
  pattern_vec <- data.frame(a=c("C","R","U"), 
                            b=c("C","R","U"), 
                            c=c("C","R","U"), 
                            d=c("C","R","U"), 
                            e=c("C","R","U")) %>% expand.grid() %>%
    unite(pattern, a:e, sep="") %>% pull(pattern)
    
  
  df_probs <- df_probs %>% group_by(chr, pos, tot_alleles) %>% arrange(pos,cat) %>%
    do(merge_codes(.$AFR, .$EUR, .$SAS, .$EAS, .$AMR, log = FALSE, patterns=pattern_vec))
  
  if (DROP_SINGLETONS) {
    write.table(df_probs, file = paste('~/Projects/rarefaction-project/data/tmp/', CHR, '_', g, '_pattern_byPosition_noSingletons.txt', sep=""),
                sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
  } else {
    write.table(df_probs, file = paste('~/Projects/rarefaction-project/data/tmp/', CHR, '_', g, '_pattern_byPosition.txt', sep=""),
                sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
  }
  
  # need to be fixed to account for lack of some patterns in the spread matrix
  df_patterns <- df_probs %>% spread(pattern, prob) %>% 
    group_by(chr) %>% 
    summarise(across(everything(), ~mean(.x))) %>% 
    gather(pattern, prob, -c(chr,pos,tot_alleles)) %>% 
    select(-c(chr,pos,tot_alleles)) %>% 
    rename(!!paste(g):=prob)
  
  if (is.null(df_all)) {
    df_all <- df_patterns
  } else {
    df_all <- inner_join(df_all, df_patterns, by="pattern")
  }
  
  # ## check if the most likely prob (sans UUUUU) matches the actual probability
  # df_rate <- df_probs %>% 
  #   group_by(pos) %>% 
  #   filter(pattern!='UUUUU') %>%
  #   filter(prob==max(prob)) %>%
  #   rename(pattern2=pattern) %>%
  #   inner_join(actual_patterns, by="pos") %>%
  #   mutate(match=(pattern==pattern2)) %>%
  #   summarise(match=sum(match)) %>%
  #   ungroup() %>%
  #   summarise(match_rate=mean(match)) %>%
  #   mutate(g=paste(g))
  # 
  # if (is.null(df_match)) {
  #   df_match <- df_rate
  # } else {
  #   df_match <- rbind(df_match, df_rate)
  # }
}

combine_pattern <- function(pattern) {
  new_pattern <- paste(str_count(pattern, 'C'), 'C', 
                       str_count(pattern, 'R'), 'R',
                       str_count(pattern, 'U'), 'U', sep='')
    return(new_pattern)
}

## g = 300 (w/ singletons) ---------------------
chr22_patterns <- read.table('~/Projects/rarefaction-project/data/tmp/300_pattern_byPosition.txt', header=T)
# set the size of the windows here
win_size_txt = '100kb'
win_size=100000
windows <- (plyr::round_any(chr22_patterns$pos - (min(chr22_patterns$pos)-1),win_size,f=ceiling) - win_size/2)/1e6
chr22_patterns$windows <- windows

## Combine Patterns - w/ UUUUU
chr22_window_pattern_vec <- chr22_patterns %>% 
  dplyr::group_by(windows) %>% 
  mutate(win_sum = sum(prob)) %>% 
  ungroup() %>% group_by(windows, pattern) %>% 
  mutate(pattern_sum = sum(prob)) %>% 
  ungroup() %>% mutate(win_prob = pattern_sum/ win_sum) %>% 
  distinct(pattern, windows, .keep_all = T) %>% 
  select(-tot_alleles, -pattern_sum, -prob, -win_sum) %>% 
  mutate(pattern=combine_pattern(pattern))

# color_pallete
color_vec <- colorRampPalette(brewer.pal(12,"Set3"))(length(unique(chr22_window_pattern_vec$pattern)))

p1 <- ggplot(chr22_window_pattern_vec %>% 
               mutate(pattern=fct_rev(pattern)), 
             aes(x=windows, y=win_prob, fill=pattern)) + 
  geom_col() + scale_fill_manual('Pattern', values=color_vec) +
  theme_pubr() + xlab('Position (Mb)') + ylab('Average Probability') +
  guides(fill=guide_legend(ncol=1))

## Combine Patterns no UUUUU
chr22_window_pattern_vec_noU <- chr22_patterns %>% 
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

p2 <- ggplot(chr22_window_pattern_vec_noU %>% 
               mutate(pattern=fct_rev(pattern)), 
             aes(x=windows, y=win_prob, fill=pattern)) + 
  geom_col() + scale_fill_manual('Pattern', values=color_vec) +
  theme_pubr() + xlab('Position (Mb)') + ylab('Average Probability') + 
  guides(fill=guide_legend(ncol=1))

## Max patterns (ungrouped)
chr22_max_patterns <- chr22_patterns %>% 
  filter(pattern!='UUUUU') %>% 
  group_by(windows) %>% 
  count(pattern) %>% 
  top_n(1)

p3 <- ggplot(chr22_max_patterns,
            aes(x=windows,
                y=fct_rev(pattern))) + 
  geom_line(aes(group=1)) +
  theme_pubr() + xlab('Position (Mb)') +
  ylab('Most Frequent Pattern')

## Max patterns (grouped)
chr22_max_patterns_grouped <- chr22_patterns %>% 
  filter(pattern!='UUUUU') %>% 
  group_by(windows) %>% 
  count(pattern) %>% 
  top_n(1) %>%
  mutate(pattern=combine_pattern(pattern))

p4 <- ggplot(chr22_max_patterns_grouped,
            aes(x=windows,
                y=pattern)) + 
  geom_line(aes(group=1)) +
  theme_pubr() + xlab('Position (Mb)') +
  ylab('Most Frequent Pattern')

p <- ggarrange(p1,p2,p4, ncol=1, common.legend = T, legend='right', align = 'v')
ggsave(p,filename = paste('~/../Downloads/chr22_patterns_', win_size_txt, '.pdf', sep=""), width = 7, height=6, units='in')


## g = 300 (no Singletons)
chr22_patterns <- read.table('~/Projects/rarefaction-project/data/tmp/22_300_pattern_byPosition_noSingletons.txt', header=T)
# set the size of the windows here
win_size_txt = '100kb'
win_size=100000
windows <- (plyr::round_any(chr22_patterns$pos - (min(chr22_patterns$pos)-1),win_size,f=ceiling) - win_size/2)/1e6
chr22_patterns$windows <- windows

## Combine Patterns - w/ UUUUU
chr22_window_pattern_vec <- chr22_patterns %>% 
  dplyr::group_by(windows) %>% 
  mutate(win_sum = sum(prob)) %>% 
  ungroup() %>% group_by(windows, pattern) %>% 
  mutate(pattern_sum = sum(prob)) %>% 
  ungroup() %>% mutate(win_prob = pattern_sum/ win_sum) %>% 
  distinct(pattern, windows, .keep_all = T) %>% 
  select(-tot_alleles, -pattern_sum, -prob, -win_sum) %>% 
  mutate(pattern=combine_pattern(pattern))

# color_pallete
color_vec <- colorRampPalette(brewer.pal(12,"Set3"))(length(unique(chr22_window_pattern_vec$pattern)))

p1 <- ggplot(chr22_window_pattern_vec %>% 
               mutate(pattern=fct_rev(pattern)), 
             aes(x=windows, y=win_prob, fill=pattern)) + 
  geom_col() + scale_fill_manual('Pattern', values=color_vec) +
  theme_pubr() + xlab('Position (Mb)') + ylab('Average Probability') +
  guides(fill=guide_legend(ncol=1))

## Combine Patterns no UUUUU
chr22_window_pattern_vec_noU <- chr22_patterns %>% 
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

p2 <- ggplot(chr22_window_pattern_vec_noU %>% 
               mutate(pattern=fct_rev(pattern)), 
             aes(x=windows, y=win_prob, fill=pattern)) + 
  geom_col() + scale_fill_manual('Pattern', values=color_vec) +
  theme_pubr() + xlab('Position (Mb)') + ylab('Average Probability') + 
  guides(fill=guide_legend(ncol=1))

## Max patterns (ungrouped)
chr22_max_patterns <- chr22_patterns %>% 
  filter(pattern!='UUUUU') %>% 
  group_by(windows) %>% 
  count(pattern) %>% 
  top_n(1)

p3 <- ggplot(chr22_max_patterns,
             aes(x=windows,
                 y=fct_rev(pattern))) + 
  geom_line(aes(group=1)) +
  theme_pubr() + xlab('Position (Mb)') +
  ylab('Most Frequent Pattern')

## Max patterns (grouped)
chr22_max_patterns_grouped <- chr22_patterns %>% 
  filter(pattern!='UUUUU') %>% 
  group_by(windows) %>% 
  count(pattern) %>% 
  top_n(1) %>%
  mutate(pattern=combine_pattern(pattern))

p4 <- ggplot(chr22_max_patterns_grouped,
             aes(x=windows,
                 y=pattern)) + 
  geom_line(aes(group=1)) +
  theme_pubr() + xlab('Position (Mb)') +
  ylab('Most Frequent Pattern')

p <- ggarrange(p1,p2,p4, ncol=1, common.legend = T, legend='right', align = 'v')
ggsave(p,filename = paste('~/../Downloads/chr22_patterns_', win_size_txt, '_noSingeltons.pdf', sep=""), width = 7, height=6, units='in')


# p1 <- ggplot(chr22_max_patterns_noUUUUU %>% 
#               group_by(windows) %>% 
#               count(pattern) %>% 
#               top_n(1) %>% 
#               ungroup() %>% 
#               mutate(group_pattern = paste(str_count(pattern, 'C'), 'C', 
#                                            str_count(pattern, 'R'), 'R', 
#                                            str_count(pattern, 'U'), 'U', sep='')),
#             aes(x=windows,
#                 y=factor(group_pattern))) + 
#   geom_line(aes(group=1))
# p1
# 
# 
# chr22_max_patterns_noUUUUU %>% mutate()
# 
# 
# 
# if (DROP_SINGLETONS) {
#   df_match_noSingletons <- df_match
# } else {
#   df_match_wSingletons <- df_match
# }
# 
# df_all2 <- df_all
# 
# ## recolor the patterns for the plot 
# keep <- df_all2 %>% gather(g, prob, -pattern) %>% 
#   filter(g==300) %>% filter(prob >= 0.01) %>%
#   pull(pattern)
# recolor_patterns <- function(pattern, keep) {
#   recolor <- ifelse(pattern %in% keep, pattern, 'Other')
#   return(data.frame(recolor=recolor))
# }
# 
# df_plot <- df_all2 %>% 
#   gather(g, prob, -pattern) %>% 
#   mutate(g=as.numeric(g)) %>%
#   do(cbind(., recolor_patterns(.$pattern, keep)))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ## plot pattern probabilities as a function of g
# mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(length(keep) + 1)
# p <- ggplot(df_plot %>% mutate(recolor = fct_reorder(recolor, prob)) %>%
#               mutate(recolor = fct_relevel(recolor, "Other")) %>% 
#               group_by(g, recolor) %>% summarise(prob=sum(prob)),
#             aes(x=g, y=prob, fill=recolor)) +
#   geom_col(lwd=0.15, color='black') + scale_fill_manual(values=mycolors) + 
#   theme_pubr(legend = 'right') + xlab("Sample Size") + ylab("Average probability of pattern") + 
#   labs(fill="Pattern") + scale_x_continuous(breaks=c(10,100,200,300), expand = c(0.005,0.005)) + scale_y_continuous(expand=c(0.005,0.005))
# 
# p
# 
# ggsave(p, filename="../../Downloads/chr22_allelePatterns_noSingletons.pdf", width=9, height=7)
# 
# 
# ## plot relative probabilities of patterns given not UUUUU -------
# 
# ## reorganize plot data
# relative_df <- df_all2 %>% 
#   gather(g, prob, -pattern) %>% 
#   mutate(g=as.numeric(g)) %>%
#   spread(g, prob) %>% 
#   filter(pattern!='UUUUU') %>% 
#   mutate(across(where(is.numeric), ~./sum(.))) %>% 
#   gather(g, prob, -pattern)
# 
# ## get actual patterns
# actual_patterns <- df_long %>% 
#   mutate(freq=minor/(minor+major)) %>% 
#   mutate(code=ifelse(freq==0, 'U', ifelse(freq <0.05, 'R', 'C'))) %>% 
#   select(-c(tot_alleles, minor, major, freq)) %>% 
#   spread(pop, code) %>% 
#   mutate(pattern=paste(AFR, EUR, SAS, EAS, AMR, sep="")) %>% 
#   select(chr, pos, pattern) %>% 
#   group_by(pattern) %>% 
#   summarise(n = n()) %>% 
#   mutate(prob=n/sum(n)) %>% select(-n) %>% 
#   mutate(g=320)
# 
# ## group only patterns >=0.005
# keep2 <- relative_df %>%
#   filter(g==300) %>%
#   filter(prob >= 0.01) %>%
#   pull(pattern)
# 
# relative_df <- relative_df %>% 
#   do(cbind(., recolor_patterns(.$pattern, keep2)))
# 
# actual_patterns <- actual_patterns %>%
#   do(cbind(., recolor_patterns(.$pattern, keep2)))
# 
# 
# levels <- actual_patterns %>% mutate(recolor=fct_reorder(recolor, prob)) %>% pull(recolor) %>% levels()
# df_plot2 <- rbind(actual_patterns, relative_df) %>%
#   mutate(g=as.numeric(g))
# 
# mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(length(levels))
# 
# p2 <- ggplot(df_plot2 %>% mutate(recolor = fct_relevel(recolor, levels)) %>% 
#                group_by(g, recolor) %>% summarise(prob=sum(prob)),
#              aes(x=g, y=prob, fill=recolor)) +
#   geom_col(lwd=0.15, color='black') + scale_fill_manual(values=mycolors) + 
#   theme_pubr(legend = 'right') + xlab("Sample Size") + ylab("Relative pattern probability") + 
#   labs(fill="Pattern") + scale_x_continuous(breaks=c(10,100,200,300), expand = c(0.005,0.005)) + scale_y_continuous(expand=c(0.005,0.005))
# 
# p2
# 
# ggsave(p2, filename="../../Downloads/chr22_allelePatterns_relative_noSingletons.pdf", width=9, height=7)
# 
# df_match <- inner_join(df_match_wSingletons %>% rename(Yes=match_rate), df_match_noSingletons %>% rename(No=match_rate)) %>% 
#   gather(`Has Singletons?`, prob, -g)
# ## plot the match frequency
# p3 <- ggplot(df_match %>% mutate(g=as.numeric(g)), aes(x=g, y=prob)) + 
#   geom_line(aes(color=`Has Singletons?`, lty=`Has Singletons?`) ,lwd=1.5) + 
#   xlab("Sample Size") + 
#   ylab("Match Frequency") +
#   theme_pubr(legend = "right") + 
#   scale_x_continuous(breaks=c(10,100,200,300), expand = c(0.005,0.005)) +
#   scale_y_continuous(limits = c(0,1), expand=c(0.005,0.005))
# p3
# 
# ggsave(p3, filename="../../Downloads/match_rate.pdf", width=4, height=4)
# 
# 
# 
# 
# x = ceiling(rnorm(5000000, mean=3000,sd=10))
# start = Sys.time()
# for (i in 1:length(x)) {choose(x[i],300)}
# end = Sys.time()
# end-start
