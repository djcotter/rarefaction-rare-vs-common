# plot_all_chrs_by_position.R
# script to take patterns per site and calculate windowed patterns across all chrs
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
              help = "input path with [] instead of the chromsome number"),
  make_option(c("--output"),
              type= "character", default = NULL,
              help = "output prefix for figure and table output"),
  make_option(c("--ext"),
              type="character", default="pdf", help="extension for figure needed for ggsave.")
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
plot_limits = c(0,250)

df_summary <- NULL
plot_list <- NULL
# loop over all chromosomes
for (i in 1:22) {
    ## read in superpop data for chr i -----
    chrom_input <- str_replace(opt$input, '\\[\\]', as.character(i))
    chr_patterns <- read.table(chrom_input, header=T)

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
    windows <- (plyr::round_any(chr_patterns$pos,
                                win_size,
                                f=ceiling))/1e6
    chr_patterns$windows <- windows

    ## Combine Patterns into ordered triples -- for averages
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

    names(myColors) <- levels

    # define plotting function
    plot_averages <- function(df, plot_colors=myColors,plot_levels=levels,chr=i) {
    if (is.null(p_limits)) {
        p_limits =  c(plyr::round_any(min(df$windows),10,f=floor),max(df$windows)+0.05)
    }
    
    x_breaks = seq(0, plyr::round_any(max(df_)), 10)
    
    p1 <- ggplot(df %>% 
                    mutate(pattern=factor(pattern,levels=plot_levels)) %>% 
                    group_by(pattern,windows) %>% 
                    select(-pos) %>% 
                    summarise(win_prob=sum(win_prob)), 
                aes(x=windows, y=win_prob, fill=pattern)) + 
        geom_col(width=win_size/1e6) + scale_fill_manual('Pattern', values=plot_colors) + 
        theme_pubr() + xlab('Position (Mb)') + ylab(paste(i)) +
        guides(fill=guide_legend(ncol=1)) + theme(legend.position='right') +
        scale_x_continuous(expand=c(0,0),limits=p_limits,breaks=x_breaks) + scale_y_continuous(expand=c(0,0),breaks=c(0,0.5,1.0)) +
        guides(fill=guide_legend(ncol = 1, override.aes = list(color="black", lwd=0.1,width=1))) +
        theme(legend.title = element_blank(), 
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.55,"line"),
            axis.title.x = element_blank(),
            axis.text = element_blank())
    
    return(p1)
    }

    p1 <- plot_averages(chr_window_pattern_vec)
    plot_list[[i]] <- p1

    ## summarise all sites into one average (as in figure 1)
    chr_patterns_summary <- chr_patterns %>% select(-pos) %>%
        group_by(pattern) %>%
        summarise(mean_prob = mean(prob)) %>%
        rename(!!paste(i) := mean_prob)
    
    if(is.null(df_summary)) {
        df_summary = chr_patterns_summary
    } else {
        df_summary = full_join(df_summary, chr_patterns_summary, by="pattern")
    }
}

for (i in seq_along(plot_list)) {
  if (i %in% 1:8) {
    plot_list[[i]] <- plot_list[[i]] + 
      scale_x_continuous(expand=c(0,0),limits=c(0,290),breaks=seq(0,290,10))
  } else {
    plot_list[[i]] <- plot_list[[i]] + 
      scale_x_continuous(expand=c(0,0),limits=c(0,140),breaks=seq(0,140,10))
  }
}

joint_legend <- get_legend(plot_list[[1]], position = "right")
group1 <- ggarrange(plotlist = plot_list[1:8], ncol = 1, align="v", legend="none")
group2 <- ggarrange(plotlist = plot_list[9:22], ncol=2, nrow=7, align="hv", legend="none")
p <- ggarrange(group1, group2, ncol=1, align="v", legend.grob = joint_legend, legend="right")
#p <- ggarrange(plotlist=plot_list, ncol=1, align="v", common.legend=T, legend="right")

scalar = 1.5
ggsave(p, filename=paste(opt$output, "_byPosition.", opt$ext, sep=""), height=11*scalar, width=8.3*scalar, units="in")

write.table(df_summary, file=paste(opt$output, '_summaryTable.txt', sep=""), sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)