# calculate_allele_patterns.R
# script to take allele counts per locus and calculate the probability of the observed pattern
# Daniel Cotter

## load packages -------
suppressPackageStartupMessages(require(tidyverse))
library(ggpubr)
library(RColorBrewer)
suppressPackageStartupMessages(require(optparse))

## parse arguments ------
option_list <- list(
  make_option(c("--input"),
              type = "character", default = NULL,
              help = "input file with all positions and all patterns probs"),
  make_option(c("--output"),
              type= "character", default = NULL,
              help = "output figure. Uses regular ggsave extensions"),
  make_option(c("--range"),
              type = "character", default = NULL,
              help = "define an integer range in which to plot the figure as START:STOP"),
  make_option(c("--rank_cutoff"),
              type="numeric", default =2,
              help = "only color ranks that ever appear at or above this cutoff in the range")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if ( is.null(opt$input) | is.null(opt$output) ) {
  print_help(opt_parser)
  stop("Must specify an input and output.", call.=FALSE)
}

## define window size and annotation -------
win_size_txt = '100kb'
win_size = 100000

## define how many of the top ranks to color ------
RANK_COLOR_CUTOFF = opt$rank_cutoff

## read in superpop data -----
chr_patterns <- read.table(opt$input, header=T)

## Define combine pattern function -------
combine_pattern <- function(pattern) {
  # define each pattern as an ordered triple
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

## define plot limits
if (is.null(opt$range)) {
  plot_limits = NULL
} else {
  plot_range = stringr::str_split(opt$range, ":", n=2)[[1]]
  plot_limits = c(as.numeric(plot_range)[1], as.numeric(plot_range[2]))
}

# set the size of the windows here
windows <- (plyr::round_any(chr_patterns$pos - (min(chr_patterns$pos)-1),win_size,f=ceiling) - win_size/2)/1e6
chr_patterns$windows <- windows

## Combine Patterns into ordered triples
chr_window_pattern_ranks <- chr_patterns %>% 
  filter(pattern!='UUUUU') %>%
  mutate(pattern=combine_pattern(pattern)) %>%
  dplyr::group_by(windows) %>% 
  mutate(win_sum = sum(prob)) %>% 
  ungroup() %>% group_by(windows, pattern) %>% 
  mutate(pattern_sum = sum(prob)) %>% 
  ungroup() %>% mutate(win_prob = pattern_sum/ win_sum) %>% 
  distinct(pattern, windows, .keep_all = T) %>% 
  select(-pattern_sum, -prob, -win_sum) %>% 
  group_by(windows) %>% 
  mutate(pattern_rank = dense_rank(desc(win_prob))) %>% 
  arrange(windows, pattern_rank)

## filter the data to fit in the plot range
if (!is.null(plot_limits)) {
  chr_window_pattern_ranks <- chr_window_pattern_ranks %>%
    filter(windows >= plot_limits[1] & windows <= plot_limits[2])
} else {
  plot_limits = c(0, max(chr_window_pattern_ranks$windows))
}

## color_pallete -----------
myColors <- c(
  "#FA7F72", "#F6BF4A",
  "#FBE3AF", "#C6B2AB",
  "#CABABC", "#CEC2CD",
  "#44AE9B", "#5CBAA9",
  "#74C6B7", "#8CD2C6",
  "#B2DE68", "#BFE47C",
  "#CCEB90", "#D9F2A3",
  "#E7F9B8", "#BB7FBD",
  "#CA96C9", "#D9ADD5",
  "#E8C3E1", "#F7DBED",
  "#F6FD91"
)

levels = c(
  "(0,0,5)", "(0,1,4)",
  "(1,0,4)", "(0,2,3)",
  "(1,1,3)", "(2,0,3)",
  "(0,3,2)", "(1,2,2)",
  "(2,1,2)", "(3,0,2)",
  "(0,4,1)", "(1,3,1)",
  "(2,2,1)", "(3,1,1)",
  "(4,0,1)", "(0,5,0)",
  "(1,4,0)", "(2,3,0)",
  "(3,2,0)", "(4,1,0)",
  "(5,0,0)"
)



names(myColors) <- levels

# drop the color for (5,0,0)
myColors <- myColors[-21]
levels <- levels[-21]

# redefine vectors for the top ranked colors
patterns_to_color = chr_window_pattern_ranks %>% filter(pattern_rank <= RANK_COLOR_CUTOFF) %>% pull(pattern) %>% unique()
all_patterns = chr_window_pattern_ranks %>% pull(pattern) %>% unique()
pattern_vec = c(patterns_to_color, all_patterns[!(all_patterns %in% patterns_to_color)])
color_vec = c(myColors[patterns_to_color], 
              rep("#D3D3D3", 20-length(patterns_to_color)))
alpha_vec = c(rep(1, length(patterns_to_color)),
              rep(0.075, 20-length(patterns_to_color))) 
lwd_vec = c(rep(0.7, length(patterns_to_color)),
            rep(0.3, 20-length(patterns_to_color))) 

## plot ranks
all_breaks=seq(20,1,-1)
my_breaks=c(20,15,10,5,1)
my_labels <- sapply(all_breaks, function(x){ ifelse(x%in%my_breaks,x,"") })

x_breaks = seq(plyr::round_any(min(chr_window_pattern_ranks$windows),10,f=floor),
               plyr::round_any(max(chr_window_pattern_ranks$windows),10, f=floor),
               10)

p_rank <- ggplot(chr_window_pattern_ranks %>%
                   mutate(pattern_rank=fct_rev(as.factor(as.numeric(pattern_rank)))), 
                 aes(x=windows, y=pattern_rank)) + 
  geom_line(aes(group=pattern, color=pattern, alpha=pattern, lwd=pattern)) +
  scale_color_manual("Pattern", breaks=pattern_vec[1:length(patterns_to_color)], values=color_vec) +
  scale_alpha_manual(breaks=pattern_vec, values=alpha_vec, guide="none") +
  scale_size_manual(breaks=pattern_vec, values=lwd_vec, guide="none") +
  xlab('Position (Mb)') +
  ylab('Rank') +
  theme_pubr(legend = "right") + 
  scale_y_discrete(breaks=all_breaks, labels=my_labels, expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0),limits=plot_limits, breaks=x_breaks) +
  guides(fill=guide_legend(ncol = 1)) +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.55,"line"))

## save the plot -----------
ggsave(p_rank,filename = opt$output, 
       width = 140, height=45, units='mm')
