libary(tidyverse)

df_wSingletons <- read.table('data/patterns/22_actualPattern_all-snps_wSingletons.txt', header=T)
df_noSingletons <- read.table('data/patterns/22_actualPattern_all-snps_noSingletons.txt',header=T)

format_data <- function(df) {
    pattern_vec <- data.frame(
        a = c("C", "R", "U"),
        b = c("C", "R", "U"),
        c = c("C", "R", "U"),
        d = c("C", "R", "U"),
        e = c("C", "R", "U")
        ) %>%
        expand.grid() %>%
        unite(pattern, a:e, sep = "") %>%
        pull(pattern)
    
    recode_pattern <- function(code_str) {
        code_map <- c("0", "1", "2")
        names(code_map) <- c("U", "R", "C")
        return(stringr::str_replace_all(code_str,code_map))
    }
    
    df2 <- left_join(data.frame(pattern=pattern_vec), df, by="pattern") %>% 
        select(pattern,n) %>%
        mutate(n = ifelse(is.na(n),0,n)) %>%
        mutate(pattern=recode_pattern(pattern))
    return(df2)
}

df_geovar_wSingletons <- format_data(df_wSingletons)
df_geovar_noSingletons <- format_data(df_noSingletons)

write.table(df_geovar_wSingletons, file="etc/22_geovar_counts_all-snps_wSingletons.txt", col.names=F, quote=F, row.names=F)
write.table(df_geovar_noSingletons, file="etc/22_geovar_counts_all-snps_noSingletons.txt", col.names=F, quote=F, row.names=F)
