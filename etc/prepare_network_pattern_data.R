library(tidyverse)


df <-  read.table('~/Projects/rarefaction-project/data/patterns/22_patterns_all-snps_noSingletons.txt', header=T)

combine_pattern <- function(pattern) {
  new_pattern <- paste(str_count(pattern, 'U'), 
                       str_count(pattern, 'R'),
                       str_count(pattern, 'C'), sep=',')
  return(new_pattern)
}

df_mod <- df %>% select(pattern, X300) %>% rename(node=pattern, prop=X300) %>% 
  mutate(node=combine_pattern(node)) %>% group_by(node) %>% summarise(prop=sum(prop))

df_mod_relative <- df_mod %>% filter(node!='5,0,0') %>% mutate(prop=prop/sum(prop))
df_mod_relative <- rbind(df_mod_relative, data.frame(node='5,0,0', prop=0))

write.table(df_mod, file = "~/Projects/rarefaction-project/etc/22_g-300_all-snps_patterns.txt", col.names=T, row.names = F, quote = F, sep='\t')
write.table(df_mod_relative, file = "~/Projects/rarefaction-project/etc/22_g-300_all-snps_relative_patterns.txt", col.names=T, row.names = F, quote = F, sep='\t')
