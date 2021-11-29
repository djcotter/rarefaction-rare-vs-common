# calculate_allele_patterns.R
# script to take allele counts per locus and calculate the probability of the observed pattern
# Daniel Cotter

## load packages -------
suppressPackageStartupMessages(require(tidyverse))
library(ggpubr)
#suppressPackageStartupMessages(require(optparse))

## temp stuff
setwd('~/Projects/rarefaction-project/')
pop_label = "superpops"
CHR = 22
g_list = seq(10,300,by = 10)
z = 0.05

## read in superpop data -----
df <- read.table(file.path('data', 'allele_counts', paste('chr', CHR, '_counts_', pop_label, '.txt', sep="")), header = TRUE)

# temporarily subset the sample size down
set.seed(1)
df_sub <- df %>% sample_n(10000, replace=FALSE)
df <- df_sub %>% tibble()

df_long <- df %>% 
  gather(pop, counts, -c(chr:major)) %>% 
  select(-c(minor:major)) %>% 
  separate(counts, into=c('minor', 'major'), sep="/") %>%
  mutate(minor=as.integer(minor), major=as.integer(major))

# remove alleles for which there is no globally "minor" allele (i.e both = 50%)
df_long <- df_long %>% filter(!is.na(minor))

## define functions to calculate the probability of (U)nobserved or (R)are -----
# (C)ommon is defined by 1 - U - R
U_prob <- function(minor, major, g) {
  g = rep(g, length(minor))
  U <- choose(major, g) / choose(minor+major, g)
  return(U)
}

R_prob <- function(minor, major, g, z) {
  g <- rep(g, length(minor))
  cutoff <- ceiling(z*(minor+major)) -1
  R <- mapply(function(N1, N2, g, k) {i=1:k; sum( choose(N1, i) * choose(N2, g-i) )}, minor, major, g, cutoff, SIMPLIFY = TRUE) / (choose(minor+major, g))
  return(R)
}

## define function to calculate probabilites of all patterns ----
merge_codes <- function(..., log=FALSE) {
  x <- data.frame(...) %>% expand.grid()
  names <- str_sub(x[[1]],1,1)
  products <- as.numeric(str_sub(x[[1]],3,-1))
  if (log) {
    for (i in 2:length(x)) {
      names <- paste(names, str_sub(x[[i]],1,1), sep="")
      products <- products + as.numeric(str_sub(x[[i]],3,-1))
    }
  } else {
    for (i in 2:length(x)) {
      names <- paste(names, str_sub(x[[i]],1,1), sep="")
      products <- products * as.numeric(str_sub(x[[i]],3,-1))
    }
  }
  return(data.frame(pattern=names, prob=products))
}

## calculate the three probabilites ----------
df_all <- NULL

for (i in 1:length(g_list)) {
  g = g_list[i]
  df_probs <- df_long %>% mutate(U = round(U_prob(minor, major, g), 6), 
                                R = round(R_prob(minor, major, g, z),6)) %>% 
    mutate(C=1-R-U) %>% 
    select(-minor, -major) %>% 
    gather(cat, prob, R, U, C) %>% 
    spread(pop, prob) 
  
  ## calculate probabilities for all patterns for each locus 
  df_probs <- df_probs %>% mutate(across(AFR:SAS, ~paste(cat, .x, sep="_"))) %>% 
    group_by(chr, pos, tot_alleles) %>% 
    do(merge_codes(.$AFR, .$AMR, .$EAS, .$EUR, .$SAS)) # can we do this without hardcoding?
  
  df_patterns <- df_probs %>% spread(pattern, prob) %>% 
    group_by(chr) %>% 
    summarise(across(everything(), ~mean(.x))) %>% 
    gather(pattern, prob, -c(chr,pos,tot_alleles)) %>% 
    select(-c(chr,pos,tot_alleles)) %>% 
    rename(!!paste(g):=prob)
  
  if (is.null(df_all)) {
    df_all <- df_patterns
    next
  }
  
  df_all <- inner_join(df_all, df_patterns, by="pattern")

}

keep <- df_all %>% filter(g==300) %>% filter(prob >= 0.01) %>% pull(pattern)
recolor_patterns <- function(pattern, keep) {
  recolor <- ifelse(pattern %in% keep, pattern, 'other')
  return(data.frame(recolor=recolor))
}

df_plot <- df_all %>% 
  gather(g, prob, -pattern) %>% 
  mutate(g=as.numeric(g)) %>%
  do(cbind(., recolor_patterns(.$pattern, keep)))


#df_plot <- df_all %>% do(cbind(., recolor_patterns(.$pattern, keep)))
  

p <- ggplot(df_plot %>% mutate(recolor = fct_relevel(recolor, "other")),
            aes(x=g, y=prob, fill=recolor)) +
  geom_col(lwd=0.5) + scale_fill_brewer(palette = 'Set3') + 
  theme_pubr(legend = 'right')

p

ggsave(p, filename="~/Downloads/chr22_allele_patterns.pdf", width=7, height=4)


