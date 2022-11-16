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
              help = "define an integer range in which to plot the figure as START:STOP")
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

# set the size of the windows here
windows <- (plyr::round_any(chr_patterns$pos - (min(chr_patterns$pos)-1),win_size,f=ceiling) - win_size/2)/1e6
chr_patterns$windows <- windows

## Combine Patterns into ordered triples
chr_window_pattern_vec <- chr_patterns %>% 
  dplyr::group_by(windows) %>% 
  mutate(win_sum = sum(prob)) %>% 
  ungroup() %>% group_by(windows, pattern) %>% 
  mutate(pattern_sum = sum(prob)) %>% 
  ungroup() %>% mutate(win_prob = pattern_sum/ win_sum) %>% 
  distinct(pattern, windows, .keep_all = T) %>% 
  select(-pattern_sum, -prob, -win_sum) %>% 
  mutate(pattern=combine_pattern(pattern))

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

names(myColors) = levels

## define plot limits
if (is.null(opt$range)) {
  plot_limits = NULL
} else {
  plot_range = stringr::str_split(opt$range, ":", n=2)[[1]]
  plot_limits = c(as.numeric(plot_range)[1], as.numeric(plot_range[2]))
}

## plot genomic positions -----
p1 <- ggplot(chr_window_pattern_vec %>% 
               mutate(pattern=factor(pattern,levels=levels)) %>% 
               group_by(pattern,windows) %>% 
               select(-pos) %>% 
               summarise(win_prob=sum(win_prob)), 
             aes(x=windows, y=win_prob, fill=pattern)) + 
  geom_col(width=win_size/1e6) + scale_fill_manual('Pattern', values=myColors) + 
  theme_pubr() + xlab('Position (Mb)') + ylab('Probability') +
  guides(fill=guide_legend(ncol=1)) + theme(legend.position='right') +
  scale_x_continuous(expand=c(0,0),limits=plot_limits) + scale_y_continuous(expand=c(0,0)) +
  guides(fill=guide_legend(ncol = 1, override.aes = list(color="black", lwd=0.1,width=1))) +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.55,"line"))

## save the file to output as a 1.5 column width
ggsave(p1,filename=opt$output, width=140,height=75, units="mm")

