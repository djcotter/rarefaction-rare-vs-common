# parse_populations.R
# script to take as input the 1kg populations and split them into individual files with lists of samples
# Daniel Cotter

## load packages -------
suppressPackageStartupMessages(require(tidyverse))

## read in data -------
df <- read.table(file.path('data', 'individual_population_codes.txt'), header = TRUE) %>%
  select(-gender)

## get pop list -------
pops <- df %>% pull(pop) %>% unique()
superpops <- df %>% pull(super_pop) %>% unique()

## create directory -----
if (!file.exists(file.path('data', 'pops'))) {
  dir.create(file.path('data', 'pops'), showWarnings = FALSE)
}


## write all samples to file -----
write.table(df %>% pull(sample),
            file=file.path('data', 'pops', 'ALL_samples.txt'),
            col.names = F, row.names = F, quote = F)

## write Superpops to file -----
for (i in 1:length(superpops)) {
  POP = superpops[i]
  pop_list = df %>% filter(super_pop == POP) %>% pull(sample)
  write.table(pop_list,
              file=file.path('data', 'pops', paste(POP, 'samples.txt', sep='_')), 
              col.names = F, row.names = F, quote = F)
}

## write Pops to file -----
for (i in 1:length(pops)) {
  POP = pops[i]
  pop_list = df %>% filter(pop == POP) %>% pull(sample)
  write.table(pop_list,
              file=file.path('data', 'pops', paste(POP, 'samples.txt', sep='_')), 
              col.names = F, row.names = F, quote = F)
}
