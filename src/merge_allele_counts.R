# merge_allele_counts.R
# script to take as input the 1kg populations and split them into individual files with lists of samples
# Daniel Cotter

## load packages -------
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(optparse))

## set seed -----
set.seed(42)

## parse arguments ------
option_list <- list(
  make_option(c("--pops"),
              type = "logical", action = "store_true", default = FALSE,
              help = "specify to merge the pops files"),
  make_option(c("--superpops"),
              type = "logical", action = "store_true", default = FALSE,
              help = "specify to merge the superpops files"),
  make_option(c("--chr"),
              type = "character", default = NULL,
              help = "chromosome to analyze and merge")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if( opt$pops && opt$superpops ) {
  print_help(opt_parser)
  stop("Can only specify --pops OR --superpops NOT BOTH.", call.=FALSE)
}

if ( !opt$pops && !opt$superpops ) {
  print_help(opt_parser)
  stop("Must specify either --pops OR --superpops.", call.=FALSE)
}

if ( is.null(opt$chr) ) {
  print_help(opt_parser)
  stop("Must specify a chromosome (e.g. --chr 22).", call.=FALSE)
}

POPS = opt$pops
SUPERPOPS = opt$superpops
CHR = paste("chr", opt$chr, sep="")

# flag to indicate whether to include or drop missing data
# i.e. data that has no observations at a loci in 1 or more populations
include_missing_data <- TRUE

## read in population data -------
pop_list <- read.table(
  file.path("data", "tmp", "individual_population_codes.txt"),
  header = TRUE) %>%
  select(-gender) %>%
  arrange(super_pop, pop)

if (POPS) {
  pop_list <- pop_list %>% pull(pop) %>% unique()
  pop_label = "pops"
} else if (SUPERPOPS) {
  pop_list <- pop_list %>% pull(super_pop) %>% unique()
  pop_label = "superpops"
}

df <- NULL
df_all <- NULL
## read in pop data ----
for (i in seq_along(pop_list)) {
  pop = pop_list[i]
  df <- read.table(file = file.path("data", "tmp", paste(pop, "_", CHR, ".frq.count", sep="")), 
                   skip=1, col.names = c("chr", "pos", "n_alleles", "tot_alleles", "allele1" , "allele2")) %>%
    tibble() %>%
    select(-n_alleles) %>% 
    gather("allele", "info", -c("chr", "pos", "tot_alleles")) %>% separate(info, into=c("allele", "count"), sep=":") %>%
    mutate(count=as.numeric(count)) %>%
    mutate(pop=pop) %>%
    mutate(freq=count/tot_alleles) %>% 
    select(-tot_alleles)
  
  if (is.null(df_all)) {
    df_all <- df
  } else {
    df_all <- rbind(df_all,df)
  }
  
}

## Handle SNPs that have data missing in one or more populations by recoding NAs in Freq column to 0
missing_data <- data.frame(pos=df_all %>% filter(is.na(freq)) %>% pull(pos) %>% unique()) %>% mutate(chr=5) %>% select(chr, pos)
write.table(missing_data,
            file=file.path("data", "allele_counts", paste(CHR, "_missing-data_", pop_label, ".txt", sep="")), 
            col.names = F, row.names = F, quote = F)

## if filtering, make a list of snps with NAs (i.e. missing data) to filter out
if (include_missing_data == FALSE) {
  missing_data <- missing_data %>% pull(pos)
} else {
  missing_data <- c()
}

# determine major and minor alleles
# na.rm if including missing data
df_all_mod <- df_all %>% 
  group_by(chr,pos) %>%
  mutate(tot_alleles=sum(count)) %>%
  ungroup() %>% 
  group_by(chr, pos, allele, tot_alleles) %>%
  summarise(avg_freq=mean(freq, na.rm=include_missing_data)) %>%
  ungroup() %>%
  mutate(class=ifelse(avg_freq<0.5,
                      'minor', ifelse(avg_freq==0.5, 
                                      'equal', 'major')))

# make a list of snps that are not variable to drop from the analysis
drop <- df_all_mod %>% filter(avg_freq==1 | avg_freq==0) %>% pull(pos) %>% unique()

## define function to randomly assign "equal" allelic types
rand_allele <- function(alleles) {
  allele1 = alleles[1]
  allele2 = alleles[2]
  test <- sample(c(0,1), 1)
  if(test){
    return(paste(allele1,allele2,sep=":"))
  } else{
    return(paste(allele2,allele1,sep=":"))
  }
}

# subset out and randomly assign one equal allele to major and one to minor
df_equal_freq <- df_all_mod %>% filter(class=="equal") %>%
  select(-avg_freq) %>%
  group_by(chr,pos,class,tot_alleles) %>% summarise(allele=rand_allele(allele)) %>%
  ungroup()

## create output directory if non existent -----
if (!file.exists(file.path("data", "allele_counts"))) {
  dir.create(file.path("data", "allele_counts"), showWarnings = FALSE)
}

write.table(df_equal_freq,
            file=file.path("data", "allele_counts", paste(CHR, "_equal-frequency-alleles_", pop_label, ".txt", sep="")), 
            col.names = F, row.names = F, quote = F)

df_all_mod <- rbind(df_all_mod %>% 
                      filter(class!="equal") %>%
                      filter(!(pos %in% missing_data)) %>%
                      select(-avg_freq) %>%
                      spread(class, allele), 
                    df_equal_freq %>%
                      select(-class) %>%
                      separate(allele,into=c("major","minor"), sep=":")) %>%
  arrange(chr,pos) %>% select(chr, pos, tot_alleles, minor, major)

## function to arrange allele counts based on previous data frame
define_minor_allele <- function(pop, a, a_count, b, b_count, minor) {
  alleles <- ifelse(a == minor, 
                    paste(a_count, b_count, sep = "/"),
                    paste(b_count, a_count, sep = "/"))
  return(alleles)
}

## loop through population files again and add counts ----
for (i in seq_along(pop_list)) {
  pop = pop_list[i]
  df <- read.table(file = file.path("data", "tmp", paste(pop, "_", CHR, ".frq.count", sep="")), 
                   skip=1, col.names = c("chr", "pos", "n_alleles", "tot_alleles", "allele1" , "allele2")) %>%
    select(-n_alleles) %>% 
    filter(!(pos %in% missing_data)) %>%
    separate(allele1, into=c("a", "a_count"), sep=":") %>% separate(allele2, into=c("b", "b_count"), sep=":")
  
  df_all_mod <- df_all_mod %>% 
    cbind(., df %>% select(a:b_count)) %>% 
    mutate(!!paste(pop) := define_minor_allele(pop, a, a_count, b, b_count, minor)) %>%
    select(-c(a, a_count, b, b_count))
}

## drop sites that are not variable
df_all_mod <- df_all_mod %>% filter(!(pos %in% drop))

## write output to file -----
write.table(df_all_mod,
            file=file.path("data", "allele_counts", paste(CHR, "_counts_", pop_label, ".txt", sep="")), 
            col.names = T, row.names = F, quote = F)