setwd('G:\\Experiments\\Paper\\FD_lifestyle_changes\\11_stat')

library(tidyverse)
library(ggplot2)
library(openxlsx)
library(ggpubr)
library(ggsci)

df <- read.xlsx(
    'G:/Experiments/Paper/FD_lifestyle_changes/ncbi_genome_table_20220812/20221107.xlsx'
)

tableS2 <- createWorkbook(creator = "ypchen", title = "wgs_strategy_influence")
# Sequencing strategy distributions
addWorksheet(tableS2,sheetName = 'wgs_strategy_distribution')
df %>%
    filter(Assembly != "GCA_024752465.1") %>%
    select(WGS_class) %>%
    table() %>%
    as.data.frame() %>%
    arrange(desc(Freq)) %>%
    mutate(proportion = round(Freq / sum(Freq) * 100, 2)) %>%
    rename(wgs_strategy = WGS_class, count = Freq) -> wgs_strategy_distribution_df
wgs_strategy_distribution_df
writeData(tableS2,sheet = 'wgs_strategy_distribution',
          x = wgs_strategy_distribution_df)
#----------------------------------------------------------------------------------


# Influence of sequencing strategies

## busco completeness

addWorksheet(tableS2,sheetName = 'completeness_nortamality')
df %>% 
    filter(Assembly != "GCA_024752465.1", WGS_class != "1st") %>%
    select(BUSCO, WGS_class) %>% group_by(WGS_class) %>% summarise(
        statistic = shapiro.test(BUSCO)$statistic,
        count = n(),
        mean = mean(BUSCO), 
        "p_value" = shapiro.test(BUSCO)$p.value,
        significance = case_when(
            shapiro.test(BUSCO)$p.value <= 0.0001 ~ "****",
            shapiro.test(BUSCO)$p.value <= 0.001 ~ "***",
            shapiro.test(BUSCO)$p.value <= 0.01 ~ "**",
            shapiro.test(BUSCO)$p.value <= 0.05 ~ "*",
            shapiro.test(BUSCO)$p.value > 0.05 ~ "ns"),
        normality = ifelse(
            shapiro.test(BUSCO)$p.value > 0.05,
            "yes",
            "no"
        )
    ) %>% 
    rename(wgs_strategy=WGS_class) -> completeness_nortamality_df

completeness_nortamality_df
writeData(tableS2,sheet = 'completeness_nortamality',
          x = completeness_nortamality_df) 

df %>%
    filter(Assembly != "GCA_024752465.1", WGS_class != "1st") %>%
    select(BUSCO, WGS_class) %>%
    ggboxplot(., x = "WGS_class", y = "BUSCO",
              color = "WGS_class", 
              add = "jitter") + 
    ylab("BUSCO completeness (%)") + xlab(NULL) +
    stat_compare_means(label.x = 1.3) +
    theme_minimal() + 
    scale_color_lancet() + theme(legend.position = 'none') -> busco_p 
busco_p
#----------------------------------------------------------

## number of contigs
addWorksheet(tableS2,sheetName = 'contig_number_nortamality')
df %>%
    filter(Assembly != "GCA_024752465.1", WGS_class != "1st") %>%
    select(Contig_scaffold_number, WGS_class) %>%
    group_by(WGS_class) %>%
    summarise(
        statistic = shapiro.test(Contig_scaffold_number)$statistic,
        count = n(),
        mean = mean(Contig_scaffold_number),
        p_value = shapiro.test(Contig_scaffold_number)$p.value,
        significance = case_when(
            shapiro.test(Contig_scaffold_number)$p.value <= 0.0001 ~ "****",
            shapiro.test(Contig_scaffold_number)$p.value <= 0.001 ~ "***",
            shapiro.test(Contig_scaffold_number)$p.value <= 0.01 ~ "**",
            shapiro.test(Contig_scaffold_number)$p.value <= 0.05 ~ "*",
            shapiro.test(Contig_scaffold_number)$p.value > 0.05 ~ "ns"),
        normality = ifelse(
            shapiro.test(Contig_scaffold_number)$p.value > 0.05,
            "yes",
            "no"
        )
    ) %>% rename(wgs_strategy=WGS_class)-> contigs_nortamality_df
contigs_nortamality_df  
writeData(tableS2,sheet = 'contig_number_nortamality',
          x = contigs_nortamality_df) 

df %>%
    filter(Assembly != "GCA_024752465.1", WGS_class != "1st") %>%
    select(Contig_scaffold_number, WGS_class) %>%
    ggboxplot(
        .,
        x = "WGS_class",
        y = "Contig_scaffold_number",
        color = "WGS_class",
        add = "jitter"
    ) +
    ylim(0, 2000) + xlab(NULL) +
    ylab("Number of Contig/Scaffold") +
    stat_compare_means(label.x = 1.3) + 
    theme_minimal() + scale_color_lancet() + 
    theme(legend.position = 'none') -> contig_number_p
contig_number_p
#-----------------------------------------------------------------------
## N50
addWorksheet(tableS2,sheetName = 'n50_normality')
df %>%
    filter(Assembly != "GCA_024752465.1", WGS_class != "1st") %>%
    select(N50, WGS_class) %>%
    group_by(WGS_class) %>%
    summarise(
        statistic = shapiro.test(N50)$statistic,
        count = n(),
        mean = mean(N50),
        p_value = shapiro.test(N50)$p.value,
        significance = case_when(
            shapiro.test(N50)$p.value <= 0.0001 ~ "****",
            shapiro.test(N50)$p.value <= 0.001 ~ "***",
            shapiro.test(N50)$p.value <= 0.01 ~ "**",
            shapiro.test(N50)$p.value <= 0.05 ~ "*",
            shapiro.test(N50)$p.value > 0.05 ~ "ns"),
        normality = ifelse(
            shapiro.test(N50)$p.value > 0.05,
            "yes",
            "no"
        )
    ) %>% rename(wgs_strategy=WGS_class) -> n50_normality_df
n50_normality_df
writeData(tableS2,sheet = 'n50_normality',
          x = n50_normality_df)

df %>%
    filter(Assembly != "GCA_024752465.1", WGS_class != "1st") %>%
    select(N50, WGS_class) %>% mutate(N50M=N50/1000000) %>%
    ggboxplot(
        .,
        x = "WGS_class",
        y = "N50M",
        color = "WGS_class",
        add = "jitter"
    ) + ylim(0,8) + 
    ylab("N50 (Mb)") + xlab(NULL) +
    stat_compare_means(label.x = 1.3, label.y = 8) + 
    theme_minimal() + scale_color_lancet() + 
    theme(legend.position = 'none') -> n50_p
n50_p

#-------------------------------------------------------------------------------

## TE size at the class level
df %>% filter(Assembly != 'GCA_024752465.1', WGS_class != "1st") %>%
    select(WGS_class, TE_size) %>% 
    group_by(WGS_class) %>%
    summarise(
        statistic = shapiro.test(TE_size)$statistic,
        count = n(),
        average = mean(TE_size),
        p_value = shapiro.test(TE_size)$p.value,
        significance = case_when(
            shapiro.test(TE_size)$p.value <= 0.0001 ~ "****",
            shapiro.test(TE_size)$p.value <= 0.001 ~ "***",
            shapiro.test(TE_size)$p.value <= 0.01 ~ "**",
            shapiro.test(TE_size)$p.value <= 0.05 ~ "*",
            shapiro.test(TE_size)$p.value > 0.05 ~ "ns"),
        normality = ifelse(
            shapiro.test(TE_size)$p.value > 0.05,
            "yes",
            "no"
        ) 
    ) %>% rename(wgs_strategy=WGS_class) -> te_size_class_df
te_size_class_df
addWorksheet(tableS2,sheetName = 'te_size_normality')
writeData(tableS2,sheet = 'te_size_normality',
          x = te_size_class_df)

df %>%
    filter(Assembly != "GCA_024752465.1", WGS_class != "1st") %>%
    select(TE_size, WGS_class) %>% mutate(TE_sizeM=TE_size/1000000) %>%
    ggboxplot(
        .,
        x = "WGS_class",
        y = "TE_sizeM",
        color = "WGS_class",
        add = "jitter"
    ) + ylim(0, 10) + xlab(NULL)+
    ylab("TE size (Mb)") + ggtitle("Sordariomycetes") +
    stat_compare_means(label.x = 1.3, label.y = 10) + 
    theme_minimal() + scale_color_lancet() + 
    theme(legend.position = 'none') -> te_size_class_p
te_size_class_p
#------------------------------------------------------------------------------------

## TE size at the family level
df %>% filter(Assembly != 'GCA_024752465.1', WGS_class != "1st") %>%
    select(WGS_class, Family) %>% 
    group_by(Family) %>% 
    table() %>% 
    as.data.frame() %>% 
    arrange(Family,desc(Freq)) %>% 
    filter(Freq>=10) %>% 
    filter(Family %in% c('Glomerellaceae', 'Hypocreaceae','Nectriaceae')) 
    
# Glomereallaceae-----------------------------------------------------------------
df %>% filter(Assembly != 'GCA_024752465.1', WGS_class != "1st") %>%
    filter(Family=='Glomerellaceae') %>%
    select(WGS_class, TE_size) %>% 
    group_by(WGS_class) %>%
    summarise(
        statistic = shapiro.test(TE_size)$statistic,
        count = n(),
        average = mean(TE_size),
        p_value = shapiro.test(TE_size)$p.value,
        significance = case_when(
            shapiro.test(TE_size)$p.value <= 0.0001 ~ "****",
            shapiro.test(TE_size)$p.value <= 0.001 ~ "***",
            shapiro.test(TE_size)$p.value <= 0.01 ~ "**",
            shapiro.test(TE_size)$p.value <= 0.05 ~ "*",
            shapiro.test(TE_size)$p.value > 0.05 ~ "ns"),
        normality = ifelse(
            shapiro.test(TE_size)$p.value > 0.05,
            "yes",
            "no"
        ) 
    ) %>% rename(wgs_strategy=WGS_class) ->te_size_glo_df 
te_size_glo_df
addWorksheet(tableS2,sheetName = 'te_size_glo_normality')
writeData(tableS2,sheet = 'te_size_glo_normality',
          x = te_size_glo_df)

df %>%
    filter(Assembly != "GCA_024752465.1", WGS_class != "1st") %>%
    filter(Family=='Glomerellaceae') %>%
    select(TE_size, WGS_class) %>% arrange(WGS_class) %>% mutate(TE_sizeM=TE_size/1000000) %>%
    ggboxplot(
        .,
        x = "WGS_class",
        y = "TE_sizeM",
        color = "WGS_class",
        add = "jitter"
    ) + ylim(0,10) +
    ylab("TE size (Mb)") + xlab(NULL) + ggtitle('Glomerellaceae') + 
    stat_compare_means(label.x = 1.3) + 
    theme_minimal() + 
    theme(legend.position = 'none') + 
    scale_color_lancet() -> te_size_glo_p
te_size_glo_p
#--------------------------------------------------------------------------------

# Nectriaceae
df %>% filter(Assembly != 'GCA_024752465.1', WGS_class != "1st") %>%
    filter(Family=='Nectriaceae') %>%
    select(WGS_class, TE_size) %>% 
    group_by(WGS_class) %>%
    summarise(
        statistic = shapiro.test(TE_size)$statistic,
        count = n(),
        average = mean(TE_size),
        p_value = shapiro.test(TE_size)$p.value,
        significance = case_when(
            shapiro.test(TE_size)$p.value <= 0.0001 ~ "****",
            shapiro.test(TE_size)$p.value <= 0.001 ~ "***",
            shapiro.test(TE_size)$p.value <= 0.01 ~ "**",
            shapiro.test(TE_size)$p.value <= 0.05 ~ "*",
            shapiro.test(TE_size)$p.value > 0.05 ~ "ns"),
        normality = ifelse(
            shapiro.test(TE_size)$p.value > 0.05,
            "yes",
            "no"
        ) 
    ) %>% rename(wgs_strategy=WGS_class) ->te_size_nec_df 
te_size_nec_df
addWorksheet(tableS2,sheetName = 'te_size_nec_normality')
writeData(tableS2,sheet = 'te_size_nec_normality',
          x = te_size_nec_df)
saveWorkbook(tableS2, file = 'tableS2_wgs_strategy.xlsx', overwrite = TRUE)

df %>%
    filter(Assembly != "GCA_024752465.1", WGS_class != "1st") %>%
    filter(Family=='Nectriaceae') %>% 
    select(TE_size, WGS_class) %>% arrange(WGS_class) %>% mutate(TE_sizeM=TE_size/1000000) %>%
    ggboxplot(
        .,
        x = "WGS_class",
        y = "TE_sizeM",
        color = "WGS_class",
        add = "jitter"
    ) +
    ggtitle("Nectriaceae") + ylab("TE size (Mb)") + xlab(NULL) +
    stat_compare_means(label.x = 1.3) + 
    theme_minimal() + 
    theme(legend.position = 'none') + 
    scale_color_lancet() -> te_size_nec_p
te_size_nec_p
#----------------------------------------------------------------------------------
ggarrange(busco_p,contig_number_p, n50_p,
          te_size_class_p, te_size_glo_p,te_size_nec_p, ncol = 3, nrow = 2, 
          labels = c('a', 'b', 'c', 'd','e','f'))

