# plot_patterns_vs_sample_size.R

# script to take pattern probabilities and plot patterns as a function
# of sample size g

# Daniel Cotter

## load packages -------
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(optparse))
library(ggpubr)
library(RColorBrewer)
library(grid)

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
opt$wSingletons <- '~/Projects/rarefaction-project/data/patterns/22_patterns_all-snps_wSingletons.txt'
opt$noSingletons <- '~/Projects/rarefaction-project/data/patterns/22_patterns_all-snps_noSingletons.txt'
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
all_colored_patterns <- rbind(df_noSingletons %>%
                                select(recolor, prob), 
                              df_wSingletons %>%
                                select(recolor, prob),
                              df_wSingletons_relative %>%
                                select(recolor, prob),
                              df_noSingletons_relative %>%
                                select(recolor, prob)
                              ) %>% 
  arrange(-prob) %>% 
  distinct(recolor,.keep_all = T)

# reorder patterns to always be in same order
levels <- all_colored_patterns %>%
  filter(recolor!='Other') %>%
  mutate(recolor=fct_reorder(recolor, prob)) %>%
  pull(recolor) %>%
  levels() %>%
  c('Other', .)

myColors <- colorRampPalette(brewer.pal(12, "Set3"))(length(levels))
names(myColors) <- levels

## define plotting functions ------
plot_patterns <- function(df_plot, colors=myColors, relative = FALSE) {
  # takes in the data frame and plots it as a function of sample size g
  plot_levels <- function(df_plot, max_g=300) {
    #max_g <- df_plot$g %>% max()
    plot_levels <- df_plot %>%
      filter(g==max_g) %>%
      group_by(g, recolor) %>%
      summarise(prob=sum(prob)) %>%
      arrange(prob) %>% filter(recolor!='Other') %>%
      pull(recolor) %>% c('Other', .)
    return(plot_levels)
  }
  pattern_levels <- plot_levels(df_plot)
  plot_colors <- function(colors, pattern_levels) {
    colors <- colors[which((names(colors) %in% pattern_levels))]
    colors <- colors[order(factor(names(colors), levels=pattern_levels))]
  }
  y_label <- if_else(relative, "Relative probability", "Average probability")
  max_g <- df_plot %>% pull(g) %>% max()
  p <- ggplot(df_plot %>%
                group_by(g, recolor) %>% 
                summarise(prob=sum(prob)) %>%
                mutate(label=ifelse(g==max_g, recolor, NA)) %>%
                mutate(recolor=fct_relevel(recolor, pattern_levels)),
              aes(x=g, y=prob, fill=recolor)) +
    geom_col(lwd=0.15, color='black') + scale_fill_manual(values=plot_colors(colors, pattern_levels)) +
    geom_text_repel(aes(label=label), position=position_stack(vjust=0.5), 
                    xlim=c(df_plot %>% pull(g) %>% max() +10, NA), 
                    size=2, direction = "y", segment.size=0.2, box.padding = 0.1,
                    force_pull=10, min.segment.length = 0.35) +
    theme_pubr(legend = 'right') + xlab("Sample size (g)") + ylab(y_label) +
    labs(fill="Pattern") + scale_x_continuous(breaks=c(10,100,200,300), expand = c(0,0), limits = c(0,max_g+60)) + scale_y_continuous(expand=c(0,0)) +
    coord_capped_cart(bottom="right")
  return(p)
}

## plot pattern vectors -------
p1 <- plot_patterns(df_wSingletons)
p2 <- plot_patterns(df_noSingletons)
p3 <- plot_patterns(df_wSingletons_relative, relative=TRUE)
p4 <- plot_patterns(df_noSingletons_relative, relative=TRUE)

legend1 <- get_legend(p1)
legend2 <- get_legend(p2)
legend3 <- get_legend(p3)
legend4 <- get_legend(p4)
legends <- ggarrange(legend1, legend3, legend2, legend4, nrow=1)
ggsave(legends, filename='~/Downloads/test_legend1.pdf', height=6)

p <- ggarrange(p1, p3, p2, p4,
               nrow = 2, ncol=2,
               widths=c(1, 1),
               labels="AUTO",
               legend = "none")
ggsave(p, filename = '~/Downloads/test_chr22_patterns.pdf', width = 190, height=170, units='mm')