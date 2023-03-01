setwd("G:\\Experiments\\Paper\\FD_lifestyle_changes\\11_stat")

library(ggplot2)
library(tidyverse)
library(ggsci)
library(readxl)
library(ggpubr)


df <- read_excel("G:\\Experiments\\Paper\\FD_lifestyle_changes\\ncbi_genome_table_20220812\\20221107.xlsx")

df_order <- read_tsv('G:\\Experiments\\Paper\\FD_lifestyle_changes\\12_ml_tree\\plot_order.tsv')
colnames(df_order) <- c('Assembly', 'plotOrder')

df_order %>% filter(Assembly != "GCA_024752465.1") -> df_order0

df %>%
  filter(Assembly != "GCA_024752465.1") %>% full_join(df_order0, by = 'Assembly') %>%
    arrange(plotOrder) %>% pull(Strain) %>% as.vector() %>% rev() -> strain_levels

df %>% 
    filter(Assembly != "GCA_024752465.1")%>% select(Strain, Lifestyle, Gene_number, Secreted_protein_number:Sucrose) -> df0




ggplot(df0, aes(x = Secreted_protein_number / Gene_number * 100, y = factor(Strain, levels = strain_levels))) +
  geom_bar(aes(fill = Lifestyle), stat = "identity", alpha=0.5) +
  theme_minimal() +
  scale_fill_lancet() + 
    xlab('Proporation of secretome (%)') + 
    ylab('Strain') + theme(legend.position = 'none')-> p1
p1

df0 %>%
  select(-Gene_number) %>%
  pivot_longer(-c(Strain, Lifestyle), names_to = "Groups", values_to = "Count") %>%
  pull(Count) %>%
  range()

print(colnames(df0))

group_level <- c(
  "Secreted_protein_number", "Effector_number", "Protease",
  "Lipase", "SSPs", "Other", "CAZy", "GH", "GT", "PL", "CE", "AA",
  "CBM", "PCWDE", "FCWDE", "Cellulose", "Hemicellulose", "Lignin",
  "Pectin", "Peptidoglycan", "Mannan", "Glucan", "Chitin", "Sucrose"
)



df0 %>%
  select(-Gene_number) %>%
  pivot_longer(-c(Strain, Lifestyle), names_to = "Groups", values_to = "Count") %>%
  ggplot(aes(x = factor(Groups, levels = group_level), y = factor(Strain, levels = strain_levels))) +
  geom_point(aes(size = Count, color = Lifestyle),
    stat = "identity",
    alpha = 0.5
  ) +
  scale_color_lancet() +
  theme_minimal() +
  scale_size(range = c(0, 10), limits = c(1, 1606)) +
  geom_text(aes(label = ifelse(Count != 0, Count, "")), size = 1) + xlab('subgroups') +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(), 
          axis.title.y = element_blank()
          )  -> p2

ggarrange(p1, p2, widths = c(1,1.5), ncol = 2)
