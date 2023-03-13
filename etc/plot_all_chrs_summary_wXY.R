suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(tidyverse))
library(ggpubr)
library(ggrepel)
library(lemon)

## parse arguments ------
option_list <- list(
  make_option(c("--input"),
              type = "character", default = NULL,
              help = "input filename"),
  make_option(c("--output"),
              type= "character", default = NULL,
              help = "output filename")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if ( is.null(opt$input) | is.null(opt$output) ) {
  print_help(opt_parser)
  stop("Must specify an input and output.", call.=FALSE)
}

# read in input file
file <- opt$input
df <- read_tsv(file)

## if sex_chromosomes
df <- df %>% rename(`23`="X", `26`="Y")
df <- df %>% pivot_longer(cols=`1`:`26`) %>%
  rename(chr=name, prob=value)

df <- df %>% pivot_longer(cols=`1`:`22`) %>%
  rename(chr=name, prob=value)

keep <- df %>% filter(chr %in% 1:22) %>% filter(prob>0.01) %>% pull(pattern) %>% unique()

recolor_patterns <- function(pattern, keep) {
  # change the pattern names to only those in the keep vector
  recolor <- ifelse(pattern %in% keep, pattern, 'Other')
  return(data.frame(recolor=recolor))
}

df <- df %>% do(cbind(., recolor_patterns(.$pattern, keep)))

levels <- df %>%
  filter(recolor!='Other') %>%
  mutate(recolor=fct_reorder(recolor, prob)) %>%
  pull(recolor) %>%
  levels() %>%
  c('Other', .)

myColors <- c(
  "#E8B063",
  "#8CD2C6",
  "#FFFFA4",
  "#D3E2D9",
  "#D2A5BE",
  "#E7F9B8",
  "#9BA7C1",
  "#C6B2AB",
  "#F6BF4A",
  "#B2DE68",
  "#EACFCD",
  "#F7DBED",
  "#CEC2CD",
  "#BB7FBD",
  "#BFCECC",
  "#FA7F72",
  "#E9E07E",
  "#96BFC4",
  "#B1A5D6",
  "#EFE7C5",
  "#F6FD91"
)
myColors <- myColors[1:19]
names(myColors) <- levels

p <- ggplot(df %>% 
              mutate(chr=as.numeric(chr)) %>%
              group_by(chr, recolor) %>% 
              summarise(prob=sum(prob)) %>%
              mutate(label=ifelse(chr==23, recolor, NA)) %>%
              mutate(recolor=fct_relevel(recolor, levels)),
            aes(x=chr, y=prob, fill=recolor)) +
  geom_col(color='black', width=1, lwd=0.05) + 
  geom_text_repel(aes(label=label), position=position_stack(vjust=0.5), 
                  xlim=c(24,28), 
                  size=1.8, direction = "y",
                  na.rm = TRUE,
                  segment.size=0.2, box.padding = 0.1,
                  force_pull=100, min.segment.length = 0.25) +
  scale_fill_manual(values=myColors) +
  theme_pubr(legend = 'none') +
  xlab("Chromosome") +
  ylab("Probability") +
  labs(fill="Pattern") + 
  scale_x_continuous(breaks=c(seq(1:23),26),
                     labels=c(seq(1:22),"X","Y"),
                     expand = c(0,0),
                     limits = c(0.5,27)) +
  scale_y_continuous(expand=c(0,0))

ggsave(p,filename=opt$output, width=7.5,height=4.5,units="in")
