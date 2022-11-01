## load packages ------
library(tidyverse)
library(ggpubr)
library(optparse)

## parse arguments ------
option_list <- list(
  make_option(
    c("--wSingletons"),
    type = "character",
    default = NULL,
    help = "path to w/ Singletons input file"
  ),
  make_option(
    c("--noSingletons"),
    type = "character",
    default = NULL,
    help = "path to no singletons input file"
  ),
  make_option(
    c("--output"),
    type = "character",
    default = NULL,
    help = "path to output file"
  )
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$wSingletons)) {
  print_help(opt_parser)
  stop("Must provide a file with patterns including singletons", call. =
         FALSE)
}

if (is.null(opt$noSingletons)) {
  print_help(opt_parser)
  stop("Must provide a file with patterns without singletons", call. = FALSE)
}

if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("Must provide an output path", call. = FALSE)
}

## create output directory if non existent ------------------
if (!file.exists(file.path("figures"))) {
  dir.create(file.path("figures"), showWarnings = FALSE)
}

## read in data ------
df_match_wSingletons <- read.table(opt$wSingletons, header=T)
df_match_noSingletons <- read.table(opt$noSingletons, header=T)

## combine data into one data frame ----
df_match <- inner_join(df_match_wSingletons %>% rename(`Singletons included`=match_rate),
                       df_match_noSingletons %>% rename(`Singletons excluded`=match_rate)) %>%
  gather(singleton, prob, -g) %>% mutate(singleton=fct_relevel(singleton, rev))
## plot the match frequency
p3 <- ggplot(df_match %>% mutate(g=as.numeric(g)), aes(x=g, y=prob)) +
  geom_line(aes(lty=singleton) ,lwd=1.5) +
  xlab("Sample size (g)") +
  ylab("Probability") +
  theme_pubr(legend = "right") +
  labs(lty="") +
  scale_x_continuous(breaks=c(10,100,200,300,400,500), expand = c(0.005,0.005)) +
  scale_y_continuous(limits = c(0,1), expand=c(0.005,0.005)) +
  theme(legend.justification = c(1,0), legend.position = c(0.9,0.1),
        legend.key.width = unit(2.5, "line"))

ggsave(p3, filename=opt$output, width=4, height=4)
