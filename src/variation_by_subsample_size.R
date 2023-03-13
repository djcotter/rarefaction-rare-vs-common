library(tidyverse)
library(ggpubr)

df <- read.table('data/patterns/22_g-500_pattern_byPosition_all-snps_noSingletons.txt', header=T)

all_positions <- df %>% pull(pos) %>% unique()

num_replicates <- 1000
subsample_size <- c(10, 100, 1000, 10000, 100000)

df_summary <- NULL
for (i in subsample_size) {
    df_all <- NULL
    for (j in 1:num_replicates) {
        subsample <- sample(all_positions, i, replace=FALSE)
        df_subsample <- df %>% filter(pos %in% subsample) %>%
                            select(-pos) %>% 
                            group_by(pattern) %>%
                            summarise(mean_prob = mean(prob)) %>%
                            arrange(pattern) %>%
                            rename(!!paste(j) := mean_prob)
        if (is.null(df_all)) {
            df_all <- df_subsample
        } else {
            df_all <- cbind(df_all, df_subsample %>% select(-pattern))
        }
    }
   df_all <- df_all %>% pivot_longer(cols=-pattern) %>% 
    group_by(pattern) %>% summarise(mean=mean(value), sd=sd(value)) %>%
    rename(!!paste("mean",i,sep="_"):=mean, !!paste("sd",i,sep="_"):=sd)
    if (is.null(df_summary)) {
        df_summary <- df_all
    } else {
        df_summary <- full_join(df_summary, df_all, by="pattern")
    }
}

df_means <- df_summary %>% unite(`10`, mean_10:sd_10) %>%
    unite(`100`, mean_100:sd_100) %>%
    unite(`1000`, mean_1000:sd_1000) %>%
    unite(`10000`, mean_10000:sd_10000) %>%
    unite(`100000`, `mean_1e+05`:`sd_1e+05`) %>%
    pivot_longer(cols=-pattern) %>%
    separate(col=value,sep="_",into=c("mean","sd")) %>%
    mutate(across(mean:sd, ~as.numeric(.))) %>%
    arrange(-mean)


keep <- df_means %>% filter(mean > 0.05) %>% pull(pattern) %>% unique()

p <- ggplot(df_means %>% filter(pattern == "RUUUU"), aes(x=name,y=mean, group=pattern, color=pattern)) + 
  geom_line() + geom_point() + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd)) +
  scale_x_log10() + theme_pubr(legend="none") + 
  labs(x="Sample Size", y="Mean Probability", title="Pattern: RUUUU")