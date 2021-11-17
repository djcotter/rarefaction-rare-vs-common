# merge_allele_counts.R
# script to take as input the 1kg populations and split them into individual files with lists of samples
# Daniel Cotter

## load packages -------
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(optparse))

## parse arguments ------
option_list = list(
  make_option(c('--pops'), type='logical', action = "store_true", default=FALSE,
              help="specify to merge the pops files"),
  make_option(c('--superpops'), type='logical', action = "store_true", default=FALSE,
              help="specify to merge the superpops files"),
  make_option(c('--chr'), type='character', default=NULL,
              help="chromosome to analyze and merge")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if( opt$pops && opt$superpops ) {
  print_help(opt_parser)
  stop("Can only specify --pops OR --superpops NOT BOTH.", call.=FALSE)
}

if( !opt$pops && !opt$superpops ) {
  print_help(opt_parser)
  stop("Must specify either --pops OR --superpops.", call.=FALSE)
}

if ( is.null(opt$chr) ) {
  print_help(opt_parser)
  stop("Must specify a chromosome (e.g. --chr 22).", call.=FALSE)
}

POPS = opt$pops
SUPERPOPS = opt$superpops
CHR = paste('chr', opt$chr, sep="")

## read in population data -------
pop_list <- read.table(file.path('data', 'tmp', 'individual_population_codes.txt'), header = TRUE) %>%
  select(-gender) %>% arrange(super_pop, pop)

if (POPS) {
  pop_list <- pop_list %>% pull(pop) %>% unique()
  pop_label = 'pops'
} else if (SUPERPOPS) {
  pop_list <- pop_list %>% pull(super_pop) %>% unique()
  pop_label = 'superpops'
}

## read in ALL data -----
df_all <- read.table(file = file.path('data', 'tmp', paste('ALL_', CHR, '.frq.count', sep="")), 
                     skip=1, col.names = c('chr', 'pos', 'n_alleles', 'tot_alleles', 'allele1' , 'allele2')) %>%
  select(-n_alleles)

## define minor allele functions -----
compare_alleles <- function(a, a_count, b, b_count) {
    alleles = ifelse(a_count == b_count, NA, 
                     ifelse(a_count<b_count, paste(a, b, sep=":"), paste(b, a, sep=":")))
    return(data.frame(alleles=alleles) %>% separate(alleles, into=c("minor", "major"), sep=":"))
}

define_minor_allele <- function(pop, a, a_count, b, b_count, minor) {
  alleles = ifelse(a == minor, paste(a_count, b_count, sep="/"), paste(b_count, a_count, sep="/"))
  # return(data.frame(alleles=alleles) %>% separate(alleles, into=c(paste(pop, "minor", sep="_"), paste(pop, "major", sep="_"))))
  return(alleles)
}

## assign globally minor alleles ------
df_all_mod <- df_all %>% separate(allele1, into=c('allele1', 'count1')) %>%
  separate(allele2, into=c('allele2', 'count2')) %>%
  mutate(count1=as.numeric(count1), count2=as.numeric(count2)) %>%
  do(cbind(., compare_alleles(.$allele1, .$count1, .$allele2, .$count2))) %>% 
  select(-allele1, -count1, -allele2, -count2)

## loop through population files ----
for (i in 1:length(pop_list)) {
  pop = pop_list[i]
  df <- read.table(file = file.path('data', 'tmp', paste(pop, "_", CHR, '.frq.count', sep="")), 
                   skip=1, col.names = c('chr', 'pos', 'n_alleles', 'tot_alleles', 'allele1' , 'allele2')) %>%
    select(-n_alleles) %>% separate(allele1, into=c('a', 'a_count'), sep=":") %>% separate(allele2, into=c('b', 'b_count'), sep=":")
  
  df_all_mod <- df_all_mod %>% 
    cbind(., df %>% select(a:b_count)) %>% 
    mutate(!!paste(pop) := define_minor_allele(pop, a, a_count, b, b_count, minor)) %>%
    select(-c(a, a_count, b, b_count))
}

## create output directory if non existent -----
if (!file.exists(file.path('data', 'allele_counts'))) {
  dir.create(file.path('data', 'allele_counts'), showWarnings = FALSE)
}

## write output to file -----
write.table(df_all_mod,
            file=file.path('data', 'allele_counts', paste(CHR, '_counts_', pop_label, '.txt', sep='')), 
            col.names = T, row.names = F, quote = F)



