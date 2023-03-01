setwd("G:\\Experiments\\Paper\\FD_lifestyle_changes\\11_stat")

library(ggplot2)
library(openxlsx)
library(tidyverse)
library(ggsci)
library(rstatix)
library(gridExtra)
library(scales)
library(ggpubr)
library(ape)
library(corrplot)
library(RColorBrewer)

df <- read.xlsx(
  "G:\\Experiments\\Paper\\FD_lifestyle_changes\\ncbi_genome_table_20220812\\20221107.xlsx"
)

# remove the outgroup and undetermined lifestyles and nematophagous fungi with only three records
df %>%
  filter(
    Assembly != "GCA_024752465.1",
    Lifestyle != "undetermined",
    Lifestyle != "nematophagous fungi"
  ) %>%
  select(Lifestyle, Gene_number, Secreted_protein_number:Sucrose) -> df0

# rename column names
df0 %>% rename(
  "Proteome" = Gene_number,
  "Secretome" = Secreted_protein_number,
  "Effector" = Effector_number
) -> df0

nrow(df0)

colnames(df0)


# get featue list, we want to present.
df0 %>%
  select(Proteome:Sucrose) %>%
  colnames() -> feature_lst
feature_lst

df0 %>%
    select(Proteome:Sucrose) -> df_corr

## crrelation analysis
addWorksheet(tableS5,sheetName = 'pearson')
colnames(df0)

p.mat <- cor.mtest(df_corr, conf.level = 0.95)$p
cor.pat = cor(df_corr)
writeData(tableS5,sheet = 'pearson',x = cor.pat, rowNames = T)

addWorksheet(tableS5,sheetName = 'pearson_significance_test')
writeData(tableS5,sheet = 'pearson_significance_test',x = p.mat, rowNames = T)


diag(cor.pat) = NA
corrplot(cor.pat, na.label = NULL,
         method = "circle", 
         diag = TRUE, 
         order = 'hclust', 
         col = brewer.pal(n = 8, name = "RdYlBu"), 
         type = "lower", tl.col = "black", 
         tl.srt = 45, 
         p.mat = p.mat, 
         insig = 'blank',
         addCoef.col = "black",
         number.cex = .5)


# for supplementary figure all functional groups--------------------------------
fig_lst <- list()
for (i in seq_along(feature_lst)) {
  if (i %in% c(1, 6, 11, 16, 21)) {
    ggplot(data = df0, aes_string(x = feature_lst[i], y = "Lifestyle", color = "Lifestyle")) +
      geom_boxplot() +
      scale_color_lancet() +
      theme_minimal() +
      theme(legend.position = "none") +
      xlab(feature_lst[i]) -> fig_lst[[i]]
  } else {
    ggplot(data = df0, aes_string(x = feature_lst[i], y = "Lifestyle", color = "Lifestyle")) +
      geom_boxplot() +
      scale_color_lancet() +
      theme_minimal() +
      theme(legend.position = "none") +
      ylab(NULL) +
      theme(axis.text.y = element_blank()) +
      xlab(feature_lst[i]) -> fig_lst[[i]]
  }
}
lay <- rbind(
  c(1, 1, 2, 3, 4, 5),
  c(6, 6, 7, 8, 9, 10),
  c(11, 11, 12, 13, 14, 15),
  c(16, 16, 17, 18, 19, 20),
  c(21, 21, 22, 23, 24, 25)
)

grid.arrange(grobs = fig_lst, layout_matrix = lay)

# for legends that need to be added to the above figure.
ggplot(data = df0, aes_string(x = feature_lst[1], y = "Lifestyle", fill = "Lifestyle")) +
  geom_boxplot() +
  scale_fill_lancet() +
  theme_minimal() +
  theme(legend.position = "top") +
  ylab(NULL) +
  theme(axis.text.y = element_blank()) +
  xlab(feature_lst[i])
#-------------------------------------------------------------------------------

## n>= 10
## n >= 10
df %>%
    filter(Assembly != "GCA_024752465.1") %>%
    group_by(Subclass) %>%
    summarise(Count = n()) %>%
    filter(Count >= 10) %>%
    pull(Subclass) %>%
    as.vector() -> n_10_subclass_lst
n_10_subclass_lst

df %>%
    filter(Assembly != "GCA_024752465.1") %>%
    group_by(Order) %>%
    summarise(Count = n()) %>%
    filter(Count >= 10) %>%
    pull(Order) %>%
    as.vector() -> n_10_order_lst
n_10_order_lst

df %>%
    filter(Assembly != "GCA_024752465.1") %>%
    group_by(Family) %>%
    summarise(Count = n()) %>%
    filter(Count >= 10) %>%
    pull(Family) %>%
    as.vector() -> n_10_family_lst
n_10_family_lst


df %>%
    filter(Assembly != "GCA_024752465.1", Lifestyle != "undetermined") %>%
    group_by(Lifestyle) %>%
    summarise(Count = n()) %>%
    filter(Count >= 10) %>%
    pull(Lifestyle) %>%
    as.vector() -> n_10_lifestyle_lst
n_10_lifestyle_lst


colnames(df0)
# difference at the taxonomic ranks

## subclass
protein_feature_lst <- feature_lst[2:25]
protein_feature_lst
tableS5 <- createWorkbook(creator = 'ypchen',title = 'protein_groups_difference')

df %>% rename("Proteome" = Gene_number,
              "Secretome" = Secreted_protein_number,
              "Effector" = Effector_number) %>%
    filter(Assembly != "GCA_024752465.1") %>%
    filter(Subclass %in% n_10_subclass_lst) %>%
    select(Subclass, Secretome:Sucrose) %>%
    group_by(Subclass) %>%
    shapiro_test(protein_feature_lst) %>%
    as.data.frame() %>%
    mutate(significance = case_when(
        p <= 0.0001 ~ "****",
        p <= 0.001 ~ "***",
        p <= 0.01 ~ "**",
        p <= 0.05 ~ "*",
        p > 0.05 ~ "ns"
    ),
    normality = ifelse(p > 0.05, "no", "yes")) %>%
    arrange(variable, Subclass) %>% 
    rename(subclass=Subclass) -> shapiro_test_subclass
shapiro_test_subclass

addWorksheet(tableS5,sheetName = 'subclass_normalitity')
writeData(tableS5,sheet = 'subclass_normalitity',x = shapiro_test_subclass)

## order
protein_feature_lst
df %>% rename("Proteome" = Gene_number,
              "Secretome" = Secreted_protein_number,
              "Effector" = Effector_number) %>%
    filter(Assembly != "GCA_024752465.1") %>%
    filter(Order %in% n_10_order_lst) %>%
    select(Order, Secretome:Sucrose) %>% 
    arrange(Order) %>%
    group_by(Order) %>%
    shapiro_test(Secretome, Effector,Protease,Lipase,SSPs,Other,CAZy,
                 GH, GT, PL, CE, AA, CBM, PCWDE, FCWDE, Cellulose,
                 Hemicellulose,Lignin,Pectin,Mannan,
                 Glucan,Chitin,Sucrose
                 ) %>%
    as.data.frame() %>%
    mutate(significance = case_when(
        p <= 0.0001 ~ "****",
        p <= 0.001 ~ "***",
        p <= 0.01 ~ "**",
        p <= 0.05 ~ "*",
        p > 0.05 ~ "ns"
    ),
    normality = ifelse(p > 0.05, "no", "yes")) %>%
    arrange(variable, Order) %>% rename(order=Order) -> shapiro_test_order
shapiro_test_order

addWorksheet(tableS5,sheetName = 'order_normalitity')
writeData(tableS5,sheet = 'order_normalitity',x = shapiro_test_order)


## family
df %>% rename("Proteome" = Gene_number,
              "Secretome" = Secreted_protein_number,
              "Effector" = Effector_number) %>%
    filter(Assembly != "GCA_024752465.1") %>%
    filter(Family %in% n_10_family_lst) %>%
    select(Family, Secretome:Sucrose) %>% 
    arrange(Family) %>%
    group_by(Family) %>%
    shapiro_test(Secretome,Effector,Protease,Lipase,SSPs,Other,CAZy,
                 GH, GT, PL,  CE, AA, PCWDE, FCWDE,Cellulose,
                 Hemicellulose,Pectin,Mannan,Glucan,Chitin
    ) %>%
    as.data.frame() %>%
    mutate(significance = case_when(
        p <= 0.0001 ~ "****",
        p <= 0.001 ~ "***",
        p <= 0.01 ~ "**",
        p <= 0.05 ~ "*",
        p > 0.05 ~ "ns"
    ),
    normality = ifelse(p > 0.05, "no", "yes")) %>%
    arrange(variable, Family) %>% rename(family=Family) -> shapiro_test_family
shapiro_test_family

addWorksheet(tableS5,sheetName = 'family_normalitity')
writeData(tableS5,sheet = 'family_normalitity',x = shapiro_test_family)

# Kruskal-Wallis Test
protein_feature_lst
## subclass
df %>% rename("Proteome" = Gene_number,
              "Secretome" = Secreted_protein_number,
              "Effector" = Effector_number) %>%
    filter(Assembly != "GCA_024752465.1") %>%
    filter(Subclass %in% n_10_subclass_lst) %>%
    select(Subclass, Secretome:Sucrose) %>%
    pivot_longer(-Subclass, names_to = "feature", values_to = "value") %>%
    group_by(feature) %>%
    kruskal_test(value ~ Subclass) %>%
    as.data.frame() %>%
    mutate(significance = case_when(p <= 0.0001 ~ '****',
                                    p <= 0.001 ~ '***',
                                    p <= 0.01 ~ '**',
                                    p <= 0.05 ~ '*',
                                    p > 0.05 ~ 'ns'),
           different = ifelse(p > 0.05, "no", "yes")) %>%
    select(-.y.)-> protein_feature_difference_subclass

protein_feature_difference_subclass

addWorksheet(tableS5,sheetName = 'subclass')
writeData(tableS5,sheet = 'subclass',x = protein_feature_difference_subclass)

## order
df %>% rename("Proteome" = Gene_number,
              "Secretome" = Secreted_protein_number,
              "Effector" = Effector_number) %>%
    filter(Assembly != "GCA_024752465.1") %>%
    filter(Order %in% n_10_order_lst) %>%
    select(Order, Secretome:Sucrose) %>%
    pivot_longer(-Order, names_to = "feature", values_to = "value") %>%
    group_by(feature) %>%
    kruskal_test(value ~ Order) %>%
    as.data.frame() %>%
    mutate(significance = case_when(p <= 0.0001 ~ '****',
                                    p <= 0.001 ~ '***',
                                    p <= 0.01 ~ '**',
                                    p <= 0.05 ~ '*',
                                    p > 0.05 ~ 'ns'),
           different = ifelse(p > 0.05, "no", "yes")) %>%
    select(-.y.)-> protein_feature_difference_order

protein_feature_difference_order

addWorksheet(tableS5,sheetName = 'order')
writeData(tableS5,sheet = 'order',x = protein_feature_difference_order)

## family
df %>% rename("Proteome" = Gene_number,
              "Secretome" = Secreted_protein_number,
              "Effector" = Effector_number) %>%
    filter(Assembly != "GCA_024752465.1") %>%
    filter(Family %in% n_10_family_lst) %>%
    select(Family, Secretome:Sucrose) %>%
    pivot_longer(-Family, names_to = "feature", values_to = "value") %>%
    group_by(feature) %>%
    kruskal_test(value ~ Family) %>%
    as.data.frame() %>%
    mutate(significance = case_when(p <= 0.0001 ~ '****',
                                    p <= 0.001 ~ '***',
                                    p <= 0.01 ~ '**',
                                    p <= 0.05 ~ '*',
                                    p > 0.05 ~ 'ns'),
           different = ifelse(p > 0.05, "no", "yes")) %>%
    select(-.y.)-> protein_feature_difference_family

protein_feature_difference_family

addWorksheet(tableS5,sheetName = 'family')
writeData(tableS5,sheet = 'family',x = protein_feature_difference_family)

## lifestyle
df %>% rename("Proteome" = Gene_number,
              "Secretome" = Secreted_protein_number,
              "Effector" = Effector_number) %>%
    filter(Assembly != "GCA_024752465.1") %>%
    filter(Lifestyle %in% n_10_lifestyle_lst) %>%
    select(Lifestyle, Secretome:Sucrose) %>%
    pivot_longer(-Lifestyle, names_to = "feature", values_to = "value") %>%
    group_by(feature) %>%
    kruskal_test(value ~ Lifestyle) %>%
    as.data.frame() %>%
    mutate(significance = case_when(p <= 0.0001 ~ '****',
                                    p <= 0.001 ~ '***',
                                    p <= 0.01 ~ '**',
                                    p <= 0.05 ~ '*',
                                    p > 0.05 ~ 'ns'),
           different = ifelse(p > 0.05, "no", "yes")) %>%
    select(-.y.)-> protein_feature_difference_lifestyle

protein_feature_difference_lifestyle

addWorksheet(tableS5,sheetName = 'lifestyle')
writeData(tableS5,sheet = 'lifestyle',x = protein_feature_difference_lifestyle)

#  pairwise difference dunntest
## subclass
df %>% rename("Proteome" = Gene_number,
              "Secretome" = Secreted_protein_number,
              "Effector" = Effector_number) %>%
    filter(Assembly != "GCA_024752465.1") %>%
    filter(Subclass %in% n_10_subclass_lst) %>%
    select(Subclass,Secretome:Sucrose) %>%
    pivot_longer(-Subclass, names_to = "feature", values_to = "value") %>%
    group_by(feature) %>%
    dunn_test(value ~ Subclass) %>%
    as.data.frame() %>%
    mutate(different = ifelse(p.adj <= 0.05, "yes", "no")) %>% 
    select(-.y.) -> protein_feature_pair_difference_subclass
protein_feature_pair_difference_subclass

df %>% rename("Proteome" = Gene_number,
              "Secretome" = Secreted_protein_number,
              "Effector" = Effector_number) %>%
    filter(Assembly != "GCA_024752465.1") %>%
    filter(Subclass %in% n_10_subclass_lst) %>%
    select(Subclass,Secretome:Sucrose) %>%
    pivot_longer(-Subclass, names_to = "feature", values_to = "value") %>%
    group_by(Subclass, feature) %>% 
    summarise(average= mean(value)) -> protein_feature_average_subclass
protein_feature_average_subclass

protein_feature_pair_difference_subclass %>% 
    left_join(protein_feature_average_subclass, 
              by = c('feature', 'group1' = 'Subclass')) %>%
    rename(group1_average=average) %>% 
    left_join(protein_feature_average_subclass, 
              by = c('feature', 'group2' = 'Subclass')) %>%
    rename(group2_average=average) -> protein_feature_pair_difference_subclass
protein_feature_pair_difference_subclass

addWorksheet(tableS5,sheetName = 'pairwise-subclass')
writeData(tableS5, sheet = 'pairwise-subclass', x = protein_feature_pair_difference_subclass)

protein_feature_pair_difference_subclass %>% 
    group_by(group1,group2) %>% 
    count(different) %>% arrange(desc(different), desc(n)) -> difference_class2_class
difference_class2_class

addWorksheet(tableS5,sheetName = 'class-class')
writeData(tableS5, sheet = 'class-class', x = difference_class2_class)

## order
df %>% rename("Proteome" = Gene_number,
              "Secretome" = Secreted_protein_number,
              "Effector" = Effector_number) %>%
    filter(Assembly != "GCA_024752465.1") %>%
    filter(Order %in% n_10_order_lst) %>%
    select(Order,Secretome:Sucrose) %>%
    pivot_longer(-Order, names_to = "feature", values_to = "value") %>%
    group_by(feature) %>%
    dunn_test(value ~ Order) %>%
    as.data.frame() %>%
    mutate(different = ifelse(p.adj <= 0.05, "yes", "no")) %>% 
    select(-.y.) -> protein_feature_pair_difference_order
protein_feature_pair_difference_order

df %>% rename("Proteome" = Gene_number,
              "Secretome" = Secreted_protein_number,
              "Effector" = Effector_number) %>%
    filter(Assembly != "GCA_024752465.1") %>%
    filter(Order %in% n_10_order_lst) %>%
    select(Order,Secretome:Sucrose) %>%
    pivot_longer(-Order, names_to = "feature", values_to = "value") %>%
    group_by(Order, feature) %>% 
    summarise(average= mean(value)) -> protein_feature_average_order
protein_feature_average_order

protein_feature_pair_difference_order %>% 
    left_join(protein_feature_average_order, 
              by = c('feature', 'group1' = 'Order')) %>%
    rename(group1_average=average) %>% 
    left_join(protein_feature_average_order, 
              by = c('feature', 'group2' = 'Order'))  %>%
    rename(group2_average=average) -> protein_feature_pair_difference_order
protein_feature_pair_difference_order






addWorksheet(tableS5,sheetName = 'pairwise-order')
writeData(tableS5, sheet = 'pairwise-order', x = protein_feature_pair_difference_order)

protein_feature_pair_difference_order %>% 
    group_by(group1,group2) %>% 
    count(different) %>% arrange(desc(different), desc(n)) -> difference_order2_order
difference_order2_order

addWorksheet(tableS5,sheetName = 'order-order')
writeData(tableS5, sheet = 'order-order', x = difference_order2_order)

## family

df %>% rename("Proteome" = Gene_number,
              "Secretome" = Secreted_protein_number,
              "Effector" = Effector_number) %>%
    filter(Assembly != "GCA_024752465.1") %>%
    filter(Family %in% n_10_family_lst) %>%
    select(Family,Secretome:Sucrose) %>%
    pivot_longer(-Family, names_to = "feature", values_to = "value") %>%
    group_by(feature) %>%
    dunn_test(value ~ Family) %>%
    as.data.frame() %>%
    mutate(different = ifelse(p.adj <= 0.05, "yes", "no")) %>% 
    select(-.y.) -> protein_feature_pair_difference_family
protein_feature_pair_difference_family

df %>% rename("Proteome" = Gene_number,
              "Secretome" = Secreted_protein_number,
              "Effector" = Effector_number) %>%
    filter(Assembly != "GCA_024752465.1") %>%
    filter(Family %in% n_10_family_lst) %>%
    select(Family,Secretome:Sucrose) %>%
    pivot_longer(-Family, names_to = "feature", values_to = "value") %>%
    group_by(Family, feature) %>% 
    summarise(average= mean(value)) -> protein_feature_average_family
protein_feature_average_family

protein_feature_pair_difference_family %>% 
    left_join(protein_feature_average_family, 
              by = c('feature', 'group1' = 'Family')) %>%
    rename(group1_average=average) %>% 
    left_join(protein_feature_average_family, 
              by = c('feature', 'group2' = 'Family'))  %>%
    rename(group2_average=average) -> protein_feature_pair_difference_family
protein_feature_pair_difference_family


addWorksheet(tableS5,sheetName = 'pairwise-family')
writeData(tableS5, sheet = 'pairwise-family', x = protein_feature_pair_difference_family)

protein_feature_pair_difference_family %>% 
    group_by(group1,group2) %>% 
    count(different) %>% arrange(desc(different), desc(n)) -> difference_family2_family
difference_family2_family

addWorksheet(tableS5,sheetName = 'family-family')
writeData(tableS5, sheet = 'family-family', x = difference_family2_family)


## lifestyle

df %>% rename("Proteome" = Gene_number,
              "Secretome" = Secreted_protein_number,
              "Effector" = Effector_number) %>%
    filter(Assembly != "GCA_024752465.1") %>%
    filter(Lifestyle %in% n_10_lifestyle_lst) %>%
    select(Lifestyle,Secretome:Sucrose) %>%
    pivot_longer(-Lifestyle, names_to = "feature", values_to = "value") %>%
    group_by(feature) %>%
    dunn_test(value ~ Lifestyle) %>%
    as.data.frame() %>%
    mutate(different = ifelse(p.adj <= 0.05, "yes", "no")) %>% 
    select(-.y.) -> protein_feature_pair_difference_lifestyle
protein_feature_pair_difference_lifestyle

df %>% rename("Proteome" = Gene_number,
              "Secretome" = Secreted_protein_number,
              "Effector" = Effector_number) %>%
    filter(Assembly != "GCA_024752465.1") %>%
    filter(Lifestyle %in% n_10_lifestyle_lst) %>%
    select(Lifestyle,Secretome:Sucrose) %>%
    pivot_longer(-Lifestyle, names_to = "feature", values_to = "value") %>%
    group_by(Lifestyle, feature) %>% 
    summarise(average= mean(value)) -> protein_feature_average_lifestyle
protein_feature_average_lifestyle

protein_feature_pair_difference_lifestyle %>% 
    left_join(protein_feature_average_lifestyle, 
              by = c('feature', 'group1' = 'Lifestyle')) %>%
    rename(group1_average=average) %>% 
    left_join(protein_feature_average_lifestyle, 
              by = c('feature', 'group2' = 'Lifestyle'))  %>%
    rename(group2_average=average) -> protein_feature_pair_difference_lifestyle

protein_feature_pair_difference_lifestyle

addWorksheet(tableS5,sheetName = 'pairwise-lifestyle')
writeData(tableS5, sheet = 'pairwise-lifestyle', x = protein_feature_pair_difference_lifestyle)

protein_feature_pair_difference_lifestyle %>% 
    group_by(group1,group2) %>% 
    count(different) %>% arrange(desc(different), desc(n)) -> difference_lifestyle2_lifestyle
difference_lifestyle2_lifestyle

addWorksheet(tableS5,sheetName = 'lifestyle-lifestyle')
writeData(tableS5, sheet = 'lifestyle-lifestyle', x = difference_lifestyle2_lifestyle)
saveWorkbook(file = 'TableS5_secretome.xlsx',tableS5, overwrite = TRUE )

#plot
##subclass
protein_feature_pair_difference_subclass %>%
    group_by(feature, different) %>%
    count() %>%
    as.data.frame() %>% 
    mutate(level='subclass') -> protein_feature_pair_difference_subclass_count_df
protein_feature_pair_difference_subclass_count_df

protein_feature_pair_difference_subclass %>%
    group_by(feature, different) %>%
    count() %>%
    as.data.frame() %>%
    ggplot(.) +
    geom_bar(aes(y = factor(feature, levels = rev(cluster_ord_lst)), 
                 x = n, 
                 fill = different), 
             position = "stack", 
             stat = "identity") +
    theme_minimal() + 
    ylab(NULL) + 
    theme(legend.position = 'none',
          axis.text.y = element_blank()) + 
    xlab(NULL) + ggtitle('Subclass')-> pair_difference_subclass_p
pair_difference_subclass_p

## order

protein_feature_pair_difference_order %>%
    group_by(feature, different) %>%
    count() %>%
    as.data.frame() %>% 
    mutate(level='order') -> protein_feature_pair_difference_order_count_df
protein_feature_pair_difference_order_count_df

protein_feature_pair_difference_order %>%
    group_by(feature, different) %>%
    count() %>%
    as.data.frame() %>%
    ggplot(.) +
    geom_bar(aes(y = factor(feature, levels = rev(cluster_ord_lst)), 
                 x = n, 
                 fill = different), 
             position = "stack", 
             stat = "identity") +
    theme_minimal() + 
    ylab(NULL) + 
    theme(legend.position = 'none',
          axis.text.y = element_blank()) + 
    xlab(NULL) + ggtitle('Order')-> pair_difference_order_p
pair_difference_order_p

## family
protein_feature_pair_difference_family %>%
    group_by(feature, different) %>%
    count() %>%
    as.data.frame() %>% 
    mutate(level='family') -> protein_feature_pair_difference_family_count_df
protein_feature_pair_difference_family_count_df

protein_feature_pair_difference_family %>%
    group_by(feature, different) %>%
    count() %>%
    as.data.frame() %>%
    ggplot(.) +
    geom_bar(aes(y = factor(feature, levels = rev(cluster_ord_lst)), 
                 x = n, 
                 fill = different), 
             position = "stack", 
             stat = "identity") +
    theme_minimal() + 
    ylab(NULL) + 
    theme(legend.position = 'none',
          axis.text.y = element_blank()) + 
    xlab(NULL) + ggtitle('Family')-> pair_difference_family_p
pair_difference_family_p


## lifestyle
protein_feature_pair_difference_lifestyle %>%
    group_by(feature, different) %>%
    count() %>%
    as.data.frame() %>% 
    mutate(level='lifestyle') -> protein_feature_pair_difference_lifestyle_count_df
protein_feature_pair_difference_lifestyle_count_df

protein_feature_pair_difference_lifestyle %>%
    group_by(feature, different) %>%
    count() %>%
    as.data.frame() %>%
    ggplot(.) +
    geom_bar(aes(y = factor(feature, levels = rev(cluster_ord_lst)), 
                 x = n, 
                 fill = different), 
             position = "stack", 
             stat = "identity") +
    theme_minimal() + 
    ylab(NULL) + 
    theme(legend.position = 'none',
          axis.text.y = element_blank()) + 
    xlab(NULL) + ggtitle('Lifestyle')-> pair_difference_lifestyle_p
pair_difference_lifestyle_p

## merge
rbind(protein_feature_pair_difference_subclass_count_df, 
      protein_feature_pair_difference_order_count_df,
      protein_feature_pair_difference_family_count_df,
      protein_feature_pair_difference_lifestyle_count_df) -> merge_df

merge_df %>% pivot_wider(names_from = c('level', 'different'), values_from = 'n') %>%
    mutate_at(c(2:9), ~replace_na(.,0)) %>% column_to_rownames('feature') -> significant_count_df
significant_count_df

addWorksheet(tableS5,sheetName = 'clustering-matrix')
writeData(tableS5, sheet = 'clustering-matrix', x = significant_count_df, rowNames = TRUE)
saveWorkbook(tableS5, "TableS5_secretome.xlsx", overwrite = TRUE)

merge_df %>% pivot_wider(names_from = c('level', 'different'), values_from = 'n') %>%
    mutate_at(c(2:9), ~replace_na(.,0)) %>% column_to_rownames('feature') %>%
    scale() %>% dist(method = 'euclidean')  %>% 
    hclust() %>% as.phylo() %>% plot() -> cluster_p
cluster_p

df %>% rename("Proteome" = Gene_number,
              "Secretome" = Secreted_protein_number,
              "Effector" = Effector_number) %>%
    filter(Assembly != "GCA_024752465.1") %>%
    filter(Lifestyle %in% n_10_lifestyle_lst) %>%
    select(Lifestyle,Secretome:Sucrose) %>% colnames()

cluster_ord_lst <- c("Peptidoglycan", "Chitin", 'CBM', 'Secretome',
                           'Other', 'Effector', 'Protease', 'Glucan', 'Cellulose',
                           'PL', 'FCWDE', 'Lignin', 'Lipase', 'GH',
                           'Mannan', 'CAZy', 'AA', 'SSPs', 'Pectin',
                           'PCWDE', 'Hemicellulose', 'CE', 'Sucrose', 'GT')

ggarrange(pair_difference_subclass_p, pair_difference_order_p, 
          pair_difference_family_p, pair_difference_lifestyle_p, 
          nrow = 1, common.legend = TRUE, widths = c(1,1,1,1)) 

# figure 7a
Top5 <- c('Cellulose', 'Hemicellulose', 'Glucan', 'FCWDE', 'Pectin')
df0
df0 %>% select(Lifestyle, all_of(Top5))-> top5df

colnames(df0)

top5df %>% pivot_longer(-c(Lifestyle), names_to = "subgroup", values_to = "count")%>% 
     group_by(Lifestyle, subgroup) %>% summarize(mean = mean(count)) %>%
    arrange(subgroup, desc(mean)) %>% filter(subgroup == 'Cellulose') %>%
    pull(Lifestyle) %>% as.vector() -> top5_order_lst

top5_order_lst
top5df$Lifestyle <- factor(top5df$Lifestyle, levels = top5_order_lst)

# figure7a
ggboxplot(top5df, x = 'Lifestyle', y = 'Cellulose', legend=FALSE, color = 'Lifestyle', 
          add = "jitter", orientation = "horizontal", size = 0.1) +
    xlab(NULL) + ylab(NULL) + ggtitle('Cellulose',) +
    geom_hline(yintercept = mean(top5df$Cellulose), linetype = 2)+ 
    stat_compare_means(method = "kruskal.test", label.y = 20, label.x = 6, size = 2)+        
    stat_compare_means(label = "p.signif", method = "wilcox.test",
                       ref.group = ".all.", hide.ns = TRUE) + 
    theme_minimal() + scale_color_lancet() + 
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) ->ps1
ps1

ggboxplot(top5df, x = 'Lifestyle', y = 'Hemicellulose', legend=FALSE, color = 'Lifestyle', 
          add = "jitter", orientation = "horizontal", size = 0.1) +
    xlab(NULL) + ylab(NULL) + ggtitle('Hemicellulose',) +
    geom_hline(yintercept = mean(top5df$Hemicellulose), linetype = 2)+ 
    stat_compare_means(method = "kruskal.test", label.y = 15, label.x = 6, size = 2)+        
    stat_compare_means(label = "p.signif", method = "wilcox.test",
                       ref.group = ".all.", hide.ns = TRUE) + 
    theme_minimal() + scale_color_lancet() + 
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
          axis.text.y = element_blank()) ->ps2
ps2

ggboxplot(top5df, x = 'Lifestyle', y = 'Glucan', legend=FALSE, color = 'Lifestyle', 
          add = "jitter" , orientation = "horizontal", size = 0.1) +
    xlab(NULL) + ylab(NULL) + ggtitle('Glucan',) +
    geom_hline(yintercept = mean(top5df$Glucan), linetype = 2)+ 
    stat_compare_means(method = "kruskal.test", label.y = 11, label.x = 6, size = 2)+        
    stat_compare_means(label = "p.signif", method = "wilcox.test",
                       ref.group = ".all.", hide.ns = TRUE) + 
    theme_minimal() + scale_color_lancet() + 
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
          axis.text.y = element_blank()) ->ps3
ps3

ggboxplot(top5df, x = 'Lifestyle', y = 'FCWDE', legend=FALSE, color = 'Lifestyle', 
          add = "jitter", orientation = "horizontal", size = 0.1) +
    xlab(NULL) + ylab(NULL) + ggtitle('FCWDE',) +
    geom_hline(yintercept = mean(top5df$FCWDE), linetype = 2)+ 
    stat_compare_means(method = "kruskal.test", label.y = 35, label.x = 6, size = 2)+        
    stat_compare_means(label = "p.signif", method = "wilcox.test",
                       ref.group = ".all.", hide.ns = TRUE) + 
    theme_minimal() + scale_color_lancet() + 
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
          axis.text.y = element_blank()) ->ps4
ps4

ggboxplot(top5df, x = 'Lifestyle', y = 'Pectin', legend=FALSE, color = 'Lifestyle', 
          add = "jitter" , orientation = "horizontal", size = 0.01, alpha = 0.5) +
    xlab(NULL) + ylab(NULL) + ggtitle('Pectin',) +
    geom_hline(yintercept = mean(top5df$Pectin), linetype = 2)+ 
    stat_compare_means(method = "kruskal.test", label.y = 15, label.x = 6, size = 2)+        
    stat_compare_means(label = "p.signif", method = "wilcox.test",
                       ref.group = ".all.", hide.ns = TRUE) + 
    theme_minimal() + scale_color_lancet(alpha = 0.5) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
          axis.text.y = element_blank()) ->ps5
ps5
ggarrange(ps1,ps2,ps3,ps4,ps5, ncol = 5, widths = c(1.6,1,1,1,1)) -> figure4b
ggsave('figure4b.pdf',figure4b, height = 2.5, width = 8.27, )

## plot 24 figures as supporting information
library(aplot)
df0 %>% select(-Proteome) -> df_24
fig24_lst <- list()
fig_tag <- df_24 %>% select(-Lifestyle) %>% colnames() %>% as.vector()
fig_tag


for ( i in seq_along(fig_tag) ) {
    tag = fig_tag[i]
    print(i)
    ggboxplot(df_24, x = 'Lifestyle', y = tag, legend=FALSE, color = 'Lifestyle', 
              add = "jitter" , orientation = "horizontal") +
        xlab(NULL) + ylab(NULL) + ggtitle(tag) +
        geom_hline(yintercept = mean(df_24[,tag]), linetype = 2)+ 
        stat_compare_means(method = "kruskal.test", 
                           label.x = 6,
                           label.y = mean(df_24[,tag]),
                           size = 2) +        
        stat_compare_means(label = "p.signif", method = "wilcox.test",
                           ref.group = ".all.", hide.ns = TRUE) + 
        theme_minimal() + scale_color_lancet(alpha = 0.5) +
        theme(legend.position = "none", 
              plot.title = element_text(size=10),
              axis.text.x = element_text(size=8),
              axis.text.y = element_blank()) -> fig24_lst[[i]]
}
    
   
ggexport(plotlist = fig24_lst, 
         filename = 'FigureS7_1.pdf', 
         width = 7.5, 
         height = 11.69, 
         nrow = 6,
         ncol = 4)




# for legend
ggboxplot(top5df, x = 'Lifestyle', y = 'Pectin', legend=FALSE, color = 'Lifestyle', 
          add = "jitter", size = 0.1) +
    xlab(NULL) + ylab(NULL) + ggtitle('Pectin',) +
    geom_hline(yintercept = mean(top5df$Pectin), linetype = 2)+ 
    stat_compare_means(method = "kruskal.test", label.y = 60)+        
    stat_compare_means(label = "p.signif", method = "wilcox.test",
                       ref.group = ".all.", hide.ns = TRUE) + 
    theme_minimal() + scale_color_lancet() + rotate_x_text(45) + 
    theme(plot.title = element_text(hjust = 0.5), legend.position = 'top') ->pslegend
pslegend

# figure4 b

figure_71_25 <- c('Cellulose', 'Hemicellulose', 'Glucan', 'Pectin',
                  'FCWDE', 'PCWDE', 'Mannan', 'Sucrose',
                  'Peptidoglycan', 'Lignin', 'Chitin', 'Effector',
                  'CE', 'AA', 'SSPs', 'PL', 
                  'CAZy', 'Secretome', 'Other', 'Lipase', 
                  'GH', 'CBM','Protease','GT',
                  'Lifestyle', 'Proteome')

df0 %>% select(all_of(figure_71_25)) -> df_71_25

figb24_lst <- list()
figb_tag <- df_71_25 %>% select(-Lifestyle, -Proteome) %>% colnames() %>% as.vector()
figb_tag
df_71_25$Lifestyle <- factor(df_71_25$Lifestyle, rev(c('entomopathogen', 
                                                   'mycoparasite', 
                                                   'saprophyte',
                                                   'endophyte', 
                                                   'plant pathogen', 
                                                   'human pathogen')))
df_71_25

df_71_25 %>% filter(Lifestyle == 'human pathogen') %>% select(Chitin) %>% 
    pull(Chitin) %>% mean()
df_71_25  %>% select(Chitin) %>% 
    pull(Chitin) %>% mean()



for ( i in seq_along(figb_tag) ) {
    tag = figb_tag[i]
    print(tag)
    ggscatter(data = df_71_25, x = 'Proteome', y = tag, 
              color = 'Lifestyle', conf.int = TRUE, size=0.5) +
        stat_cor(aes(color = Lifestyle), size = 1.75) + 
        xlab(NULL) + 
        ylab(tag) +
        geom_smooth(aes_string(x = "Proteome" , y= tag), method = 'lm') + 
        theme_minimal() + 
        scale_color_lancet() + 
        ggtitle(tag) + 
        theme(legend.position = 'none',
              plot.title = element_text(size=10),
              axis.text.x = element_text(size=8),
              axis.title.x = element_text(size=10),
              axis.title.y = element_text(size=10),
        ) -> figb24_lst[[i]]
}
figb24_lst[[i]]
ggexport(plotlist = figb24_lst, 
         filename = 'FigureS7_2.pdf', 
         width = 8.27, 
         height = 11.69, 
         nrow = 6,
         ncol = 4)

help(geom_smooth)

# check legend
ggscatter(data = df_71_25, x = 'Proteome', y = 'Cellulose', 
          color = 'Lifestyle', conf.int = TRUE, size=0.5) +
    stat_cor(aes(color = Lifestyle), size = 1.75) + 
    xlab(NULL) + 
    ylab('Cellulose') +
    geom_smooth(aes_string(x = "Proteome" , y= 'Cellulose'), method = 'lm') + 
    theme_minimal() + 
    scale_color_lancet() + 
    ggtitle(tag) + 
    theme(legend.position = 'top',
          plot.title = element_text(size=10),
          axis.text.x = element_text(size=8),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
    )

help(stat_cor)





ggscatter(data = df7, y = 'Secretome', x = 'Hemicellulose', 
          color = 'Lifestyle', palette = "jco", conf.int = TRUE, size=0.5) +
    stat_cor(aes(color = Lifestyle), label.x = 0, size = 1.5) + 
    xlab(NULL) + 
    ylab(NULL) +
    geom_smooth(aes(x=Hemicellulose,y=Secretome), method = 'lm') + 
    theme_minimal() + 
    scale_color_lancet() + 
    ggtitle('Hemicellulose') + 
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = 'none', axis.text.y = element_blank()) ->figure4c2
figure4c2

ggscatter(data = df7, y = 'Secretome', x = 'Glucan', 
          color = 'Lifestyle', palette = "jco", conf.int = TRUE, size=0.5) +
    stat_cor(aes(color = Lifestyle), label.x = 0, size = 1.5) + 
    xlab(NULL) + 
    ylab(NULL) +
    geom_smooth(aes(x=Glucan,y=Secretome), method = 'lm') + 
    theme_minimal() + 
    scale_color_lancet() + 
    ggtitle('Glucan') + 
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = 'none', axis.text.y = element_blank()) ->figure4c3
figure4c3

ggscatter(data = df7, y = 'Secretome', x = 'FCWDE', 
          color = 'Lifestyle', palette = "jco", conf.int = TRUE, size=0.5) +
    stat_cor(aes(color = Lifestyle), label.x = 8, size = 1.5) + 
    xlab(NULL) + 
    ylab(NULL) +
    geom_smooth(aes(x=FCWDE,y=Secretome), method = 'lm') + 
    theme_minimal() + 
    scale_color_lancet() + 
    ggtitle('FCWDE') + 
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = 'none', axis.text.y = element_blank()) ->figure4c4
figure4c4

ggscatter(data = df7, y = 'Secretome', x = 'Pectin', 
          color = 'Lifestyle', palette = "jco", conf.int = TRUE, size=0.5) +
    stat_cor(aes(color = Lifestyle), label.x = 0, size = 1.5) + 
    xlab(NULL) + 
    ylab(NULL) +
    geom_smooth(aes(x=Pectin,y=Secretome), method = 'lm') + 
    theme_minimal() + 
    scale_color_lancet() + 
    ggtitle('Pectin') + 
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = 'none', axis.text.y = element_blank()) ->figure4c5
figure4c5

ggarrange(figure4c1,figure4c2,figure4c3,figure4c4,figure4c5, 
          ncol = 5, widths = c(1.3,1,1,1,1)) -> figure4c

ggsave('figure4c.pdf',figure4c, height = 2.5, width = 8.27, )
# for 4c legend
ggscatter(data = df7, y = 'Secretome', x = 'Pectin', 
          color = 'Lifestyle', palette = "jco", conf.int = TRUE, size=0.5) +
    stat_cor(aes(color = Lifestyle), label.x = 0) + 
    xlab(NULL) + 
    ylab(NULL) +
    geom_smooth(aes(x=Pectin,y=Secretome), method = 'lm') + 
    theme_minimal() + 
    scale_color_lancet() + 
    ggtitle('Pectin') + 
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = 'top', axis.text.y = element_blank()) ->figure4clegend
figure4clegend


# figure4 d
correaltion_7 <- c('Cellulose', 'Hemicellulose', 'Glucan', 'FCWDE', 
                   'Pectin', 'Proteome', 'Secretome', 'Lifestyle')
ggscatter(data = df7, y = 'Proteome', x = 'Cellulose', 
          color = 'Lifestyle', palette = "jco", conf.int = TRUE, size=0.5) +
    stat_cor(aes(color = Lifestyle), label.x = 0, size = 1.5) + 
    xlab(NULL) + 
    ylab('Proteome') +
    geom_smooth(aes(y=Proteome, x=Cellulose), method = 'lm') + 
    theme_minimal() + 
    scale_color_lancet() + 
    ggtitle('Cellulose') + 
    theme(plot.title = element_text(hjust = 0.5), legend.position = 'none') -> figure4d1
figure4d1
ggscatter(data = df7, y = 'Proteome', x = 'Hemicellulose', 
          color = 'Lifestyle', palette = "jco", conf.int = TRUE, size=0.5) +
    stat_cor(aes(color = Lifestyle), label.x = 0, size = 1.5) + 
    xlab(NULL) + 
    ylab(NULL) +
    geom_smooth(aes(y=Proteome, x=Hemicellulose), method = 'lm') + 
    theme_minimal() + 
    scale_color_lancet() + 
    ggtitle('Hemicellulose') + 
    theme(plot.title = element_text(hjust = 0.5), legend.position = 'none',
          axis.text.y = element_blank()) -> figure4d2
figure4d2
ggscatter(data = df7, y = 'Proteome', x = 'Glucan', 
          color = 'Lifestyle', palette = "jco", conf.int = TRUE, size=0.5) +
    stat_cor(aes(color = Lifestyle), label.x = 0, size = 1.5) + 
    xlab(NULL) + 
    ylab(NULL) +
    geom_smooth(aes(y=Proteome, x=Glucan), method = 'lm') + 
    theme_minimal() + 
    scale_color_lancet() + 
    ggtitle('Glucan') + 
    theme(plot.title = element_text(hjust = 0.5), legend.position = 'none',
          axis.text.y = element_blank()) -> figure4d3
figure4d3
ggscatter(data = df7, y = 'Proteome', x = 'FCWDE', 
          color = 'Lifestyle', palette = "jco", conf.int = TRUE, size=0.5) +
    stat_cor(aes(color = Lifestyle), label.x = 0, size = 1.5) + 
    xlab(NULL) + 
    ylab(NULL) +
    geom_smooth(aes(y=Proteome, x=FCWDE), method = 'lm') + 
    theme_minimal() + 
    scale_color_lancet() + 
    ggtitle('FCWDE') + 
    theme(plot.title = element_text(hjust = 0.5), legend.position = 'none',
          axis.text.y = element_blank()) -> figure4d4
figure4d4
ggscatter(data = df7, y = 'Proteome', x = 'Pectin', 
          color = 'Lifestyle', palette = "jco", conf.int = TRUE, size=0.5) +
    stat_cor(aes(color = Lifestyle), label.x = 0, size = 1.5) + 
    xlab(NULL) + 
    ylab(NULL) +
    geom_smooth(aes(y=Proteome, x=Pectin), method = 'lm') + 
    theme_minimal() + 
    scale_color_lancet() + 
    ggtitle('Pectin') + 
    theme(plot.title = element_text(hjust = 0.5), legend.position = 'none',
          axis.text.y = element_blank()) -> figure4d5
figure4d5
ggarrange(figure4d1,figure4d2,figure4d3,figure4d4,figure4d5, 
          ncol = 5, widths = c(1.3,1,1,1,1)) -> figure4d
figure4d
ggsave('figure4d.pdf',figure4d, height = 2.5, width = 8.27, )

# for figure4d legend
ggscatter(data = df7, y = 'Proteome', x = 'Pectin', 
          color = 'Lifestyle', palette = "jco", conf.int = TRUE, size=0.5) +
    stat_cor(aes(color = Lifestyle), label.x = 0, size = 1.5) + 
    xlab(NULL) + 
    ylab(NULL) +
    geom_smooth(aes(y=Proteome, x=Pectin), method = 'lm') + 
    theme_minimal() + 
    scale_color_lancet() + 
    ggtitle('Pectin') + 
    theme(plot.title = element_text(hjust = 0.5), legend.position = 'top',
          axis.text.y = element_blank()) -> figure4d_legend
figure4d_legend






