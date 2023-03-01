getwd()
setwd("G:\\Experiments\\Paper\\FD_lifestyle_changes\\11_stat")
dir()

library(ggpubr)
library(treedataverse)
library(tidyverse)
library(ggsci)
library(readxl)
library(ggplot2)

df <- read_excel(
    "G:\\Experiments\\Paper\\FD_lifestyle_changes\\11_stat\\TableS_TEs_abundance.xlsx")

colnames(df)
# check range
df %>% pivot_longer(-c(Assembly:Family, Genome_size, TE_content), names_to = 'subgroups', values_to = 'base_pair') %>%
    mutate(proporation=base_pair / Genome_size * 100) %>% pull(proporation) %>% range()

# plot order
df %>% arrange(Tree_order) %>% pull(Strain) %>% as.vector() %>% rev() -> plot_order

# geom_text > 0.01
df %>% pivot_longer(-c(Assembly:Family, Genome_size, TE_content), names_to = 'subgroups', values_to = 'base_pair') %>%
    mutate(proporation=base_pair / Genome_size * 100) %>% filter(proporation >=0.01) -> text_subset

df %>% pivot_longer(-c(Assembly:Family, Genome_size, TE_content), names_to = 'subgroups', values_to = 'base_pair') %>%
    mutate(proporation=base_pair / Genome_size * 100) %>%
    ggplot(., aes(y=factor(Strain, levels = plot_order), x=subgroups)) + 
    geom_point(aes(size = proporation, color=Lifestyle), alpha=00.75) + 
    scale_color_lancet() + scale_size(range =c(3,18), limits = c(0.01,60)) + 
    theme_minimal() + 
    geom_text(data =text_subset, aes(label=round(base_pair / Genome_size * 100, 2)), color="white", size = 1.5) +
    scale_x_discrete(position = "top") -> bubble_p

df %>% ggplot(.) +
    geom_bar(aes(y=factor(Strain, levels = plot_order), x = TE_content, fill = Lifestyle), stat = "identity", alpha= 0.85) + 
    scale_fill_lancet() + 
    theme_minimal() + scale_x_continuous(position = 'top') + 
    theme(axis.text.y = element_blank(), legend.position = 'none') + ylab(NULL) -> bar_p

ggarrange(bubble_p,bar_p, widths = c(7,1))
    
    



