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
g = 500
z = 0.05
DROP_SINGLETONS = TRUE

## read in superpop data -----
df <- read.table('~/Projects/rarefaction-project/data/patterns/6_g-500_pattern_byPosition_all-snps_noSingletons.txt', header=T)

combine_pattern <- function(pattern) {
  new_pattern <- paste(
    '(',
    str_count(pattern, 'U'),
    ',',
    str_count(pattern, 'R'),
    ',',
    str_count(pattern, 'C'),
    ')',
    sep = ''
  )
  return(new_pattern)
}

## g = 300 (w/ singletons) ---------------------
chr_patterns <- df
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
  select(-pattern_sum, -prob, -win_sum) %>% 
  mutate(pattern=combine_pattern(pattern))

# color_pallete
levels <- chr_window_pattern_vec$pattern %>% unique()

# myColors <- c(
#   "#E8B063",
#   "#8CD2C6",
#   "#FFFFA4",
#   "#D3E2D9",
#   "#D2A5BE",
#   "#E7F9B8",
#   "#9BA7C1",
#   "#C6B2AB",
#   "#F6BF4A",
#   "#B2DE68",
#   "#EACFCD",
#   "#F7DBED",
#   "#CEC2CD",
#   "#BB7FBD",
#   "#BFCECC",
#   "#FA7F72",
#   "#E9E07E",
#   "#96BFC4",
#   "#B1A5D6",
#   "#EFE7C5",
#   "#F6FD91"
# )

myColors<- c(
  "#FA7F72",
  "#F6BF4A",
  "#FBE3AF",
  "#C6B2AB",
  "#CABABC",
  "#CEC2CD",
  "#44AE9B",
  "#5CBAA9",
  "#74C6B7",
  "#8CD2C6",
  "#B2DE68",
  "#BFE47C",
  "#CCEB90",
  "#D9F2A3",
  "#E7F9B8",
  "#BB7FBD",
  "#CA96C9",
  "#D9ADD5",
  "#E8C3E1",
  "#F7DBED",
  "#F6FD91"
)

myColors %>% scales::show_col()

names(myColors) = levels

p1 <- ggplot(chr_window_pattern_vec %>% 
               mutate(pattern=factor(pattern,levels=levels)) %>% 
               group_by(pattern,windows) %>% 
               select(-pos) %>% 
               summarise(win_prob=sum(win_prob)), 
             aes(x=windows, y=win_prob, fill=pattern)) + 
  geom_col(width=win_size/1e6) + scale_fill_manual('Pattern', values=myColors) + 
  theme_pubr() + xlab('Position (Mb)') + ylab('Probability') +
  guides(fill=guide_legend(ncol=1)) + theme(legend.position='right') +
  scale_x_continuous(expand=c(0,0),limits=c(28,35)) + scale_y_continuous(expand=c(0,0)) +
  guides(fill=guide_legend(ncol = 1, override.aes = list(color="black", lwd=0.1,width=1))) +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.55,"line"))
ggsave(p1,filename='~/../Downloads/test_chr6_byPosition.pdf', width=140,height=75, units="mm")

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

