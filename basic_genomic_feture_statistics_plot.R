setwd("G:\\Experiments\\Paper\\FD_lifestyle_changes\\11_stat")

library(tidyverse)
library(openxlsx)
library(ggplot2)
library(ggcorrplot)
library(broom)
library(rstatix)
library(openxlsx)
library(ggsci)
library(RColorBrewer)
library(corrplot)
library(ape)
library(scales)
library(ggpubr)

df <- read.xlsx(
  "G:\\Experiments\\Paper\\FD_lifestyle_changes\\ncbi_genome_table_20220812\\20221107.xlsx"
)

colnames(df)

df %>%
  filter(Assembly != "GCA_024752465.1") %>%
  select(Genome_size:tRNA_length_mean, -c(Contig_scaffold_number, N50, BUSCO, GO_number, Protein_number)) %>%
  select(
    Gene_number, Gene_length_min, Gene_length_mean, Gene_length_max,
    Exon_number, Exon_length_min, Exon_length_mean, Exon_length_max,
    Intron_number, Intron_length_min, Intron_length_mean, Intron_length_max,
    Intergenic_length_min, Intergenic_length_mean, Intergenic_length_max,
    tRNA_number, tRNA_length_min, tRNA_length_mean, tRNA_length_max,
    GC_content_TE, GC_content_without_TE, GC_content,
    TE_size, Genome_size_without_TE, Genome_size,
  ) -> df0
colnames(df0)

# figure3 correlation---------------------------------------------------------------


tableS3 <- createWorkbook(creator = "ypchen", title = "feature_correlation_analysis")
addWorksheet(tableS3,sheetName = 'pearson')

p.mat <- cor.mtest(df0, conf.level = 0.95)$p
cor.pat = cor(df0)
writeData(tableS3,sheet = 'pearson',x = cor.pat, rowNames = T)

addWorksheet(tableS3,sheetName = 'significance_test')
writeData(tableS3,sheet = 'significance_test',x = p.mat, rowNames = T)
saveWorkbook(tableS3, file = 'tableS3_basic_genomic_correlation_analysis.xlsx', overwrite = TRUE)

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
#------------------------------------------------------------------------------

# compare difference at taxonomic levels and different lifestyles
## n >= 10
df %>%
  filter(Assembly != "GCA_024752465.1") %>%
  group_by(Subclass) %>%
  summarise(Count = n()) %>%
  filter(Count >= 10) %>%
  pull(Subclass) %>%
  as.vector() -> n_10_subclass_lst

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

df %>%
  filter(Assembly != "GCA_024752465.1", Lifestyle != "undetermined") %>%
  group_by(Lifestyle) %>%
  summarise(Count = n()) %>%
  filter(Count >= 10) %>%
  pull(Lifestyle) %>%
  as.vector() -> n_10_lifestyle_lst
n_10_lifestyle_lst

## normality test
# out text table
tableS4 <- createWorkbook(creator = 'ypchen',title = 'genomic_feature_difference')
df %>%
  filter(Assembly != "GCA_024752465.1") %>%
  filter(Subclass %in% n_10_subclass_lst) %>%
  select(Subclass, Genome_size:tRNA_length_mean, -c(Contig_scaffold_number, N50, BUSCO, Protein_number, GO_number)) %>%
  group_by(Subclass) %>%
  shapiro_test(
    Genome_size, Genome_size_without_TE, Genome_size_without_TE,
    GC_content, GC_content_without_TE, GC_content_TE,
    tRNA_length_max, tRNA_length_mean, tRNA_length_min, tRNA_number,
    Intergenic_length_max, Intergenic_length_mean, Intergenic_length_min,
    Intron_length_max, Intron_length_mean, Intron_length_min, Intron_number,
    Exon_length_max, Exon_length_mean, Exon_number,
    Gene_length_max, Gene_length_mean, Gene_length_min, Gene_number
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
  arrange(variable, Subclass) %>% rename(subclass=Subclass) -> shapiro_test_subclass
shapiro_test_subclass

addWorksheet(tableS4,sheetName = 'subclass_normalitity')
writeData(tableS4,sheet = 'subclass_normalitity',x = shapiro_test_subclass)

df %>%
  filter(Assembly != "GCA_024752465.1") %>%
  filter(Order %in% n_10_order_lst) %>%
  select(Order, Genome_size:tRNA_length_mean, -c(Contig_scaffold_number, N50, BUSCO, Protein_number, GO_number)) %>%
  group_by(Order) %>%
  shapiro_test(
    Genome_size, Genome_size_without_TE, Genome_size_without_TE,
    GC_content, GC_content_without_TE, GC_content_TE,
    tRNA_length_max, tRNA_length_mean, tRNA_length_min, tRNA_number,
    Intergenic_length_max, Intergenic_length_mean, Intergenic_length_min,
    Intron_length_max, Intron_length_mean, Intron_length_min, Intron_number,
    Exon_length_max, Exon_length_mean, Exon_number,
    Gene_length_max, Gene_length_mean, Gene_length_min, Gene_number
  ) %>%
  as.data.frame() %>%
  mutate(significance = case_when(
      p <= 0.0001 ~ "****",
      p <= 0.001 ~ "***",
      p <= 0.01 ~ "**",
      p <= 0.05 ~ "*",
      p > 0.05 ~ "ns"
  ),normality = ifelse(p > 0.05, "yes", "no")) %>%
  arrange(variable, Order) %>% rename(order=Order)-> shapiro_test_order
shapiro_test_order
addWorksheet(tableS4,sheetName = 'order_normalitity')
writeData(tableS4,sheet = 'order_normalitity',x = shapiro_test_order)

df %>%
  filter(Assembly != "GCA_024752465.1") %>%
  filter(Family %in% n_10_family_lst) %>%
  select(Family, Genome_size:tRNA_length_mean, -c(Contig_scaffold_number, N50, BUSCO, Protein_number, GO_number)) %>%
  group_by(Family) %>%
  shapiro_test(
    Genome_size, Genome_size_without_TE, Genome_size_without_TE,
    GC_content, GC_content_without_TE, GC_content_TE,
    tRNA_length_max, tRNA_length_mean, tRNA_length_min, tRNA_number,
    Intergenic_length_max, Intergenic_length_mean, Intergenic_length_min,
    Intron_length_max, Intron_length_mean, Intron_length_min, Intron_number,
    Exon_length_max, Exon_length_mean, Exon_number,
    Gene_length_max, Gene_length_mean, Gene_length_min, Gene_number
  ) %>%
  as.data.frame() %>%
  mutate(significance = case_when(
      p <= 0.0001 ~ "****",
      p <= 0.001 ~ "***",
      p <= 0.01 ~ "**",
      p <= 0.05 ~ "*",
      p > 0.05 ~ "ns"
  ),normality = ifelse(p < 0.05, "yes", "no")) %>%
  arrange(variable, Family) %>% rename(family=Family)-> shapiro_test_family
shapiro_test_family

addWorksheet(tableS4,sheetName = 'family_normalitity')
writeData(tableS4,sheet = 'family_normalitity',x = shapiro_test_family)

df %>%
  filter(Assembly != "GCA_024752465.1", Lifestyle != "undetermined") %>%
  filter(Lifestyle %in% n_10_lifestyle_lst) %>%
  select(Lifestyle, Genome_size:tRNA_length_mean, -c(Contig_scaffold_number, N50, BUSCO, Protein_number, GO_number)) %>%
  group_by(Lifestyle) %>%
  shapiro_test(
    Genome_size, Genome_size_without_TE, Genome_size_without_TE,
    GC_content, GC_content_without_TE, GC_content_TE,
    tRNA_length_max, tRNA_length_mean, tRNA_length_min, tRNA_number,
    Intergenic_length_max, Intergenic_length_mean, Intergenic_length_min,
    Intron_length_max, Intron_length_mean, Intron_length_min, Intron_number,
    Exon_length_max, Exon_length_mean, Exon_number,
    Gene_length_max, Gene_length_mean, Gene_length_min, Gene_number
  ) %>%
  as.data.frame() %>%
  mutate(significance = case_when(
      p <= 0.0001 ~ "****",
      p <= 0.001 ~ "***",
      p <= 0.01 ~ "**",
      p <= 0.05 ~ "*",
      p > 0.05 ~ "ns"
  ),normality = ifelse(p < 0.05, "yes", "no")) %>%
  arrange(variable, Lifestyle) %>% rename(lifestyle=Lifestyle)-> shapiro_test_lifestyle
shapiro_test_lifestyle

addWorksheet(tableS4,sheetName = 'lifestyle_normalitity')
writeData(tableS4,sheet = 'lifestyle_normalitity',x = shapiro_test_lifestyle)

# Kruskal-Wallis Test
## subclass
feature_lst <- c(
  "Genome_size", "Genome_size_without_TE", "Genome_size_without_TE",
  "GC_content", "GC_content_without_TE", "GC_content_TE",
  "tRNA_length_max", "tRNA_length_mean", "tRNA_length_min", "tRNA_number",
  "Intergenic_length_max", "Intergenic_length_mean", "Intergenic_length_min",
  "Intron_length_max", "Intron_length_mean", "Intron_length_min", "Intron_number",
  "Exon_length_max", "Exon_length_mean", "Exon_number",
  "Gene_length_max", "Gene_length_mean", "Gene_length_min", "Gene_number"
)

df %>%
  filter(Assembly != "GCA_024752465.1") %>%
  filter(Subclass %in% n_10_subclass_lst) %>%
  select(Subclass, Genome_size:tRNA_length_mean, -c(Contig_scaffold_number, N50, BUSCO, Protein_number, GO_number)) %>%
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
    select(-.y.)-> geomic_feature_difference_subclass

geomic_feature_difference_subclass

addWorksheet(tableS4,sheetName = 'subclass')
writeData(tableS4,sheet = 'subclass',x = geomic_feature_difference_subclass)

## order
df %>%
  filter(Assembly != "GCA_024752465.1") %>%
  filter(Order %in% n_10_order_lst) %>%
  select(Order, Genome_size:tRNA_length_mean, -c(Contig_scaffold_number, N50, BUSCO, Protein_number, GO_number)) %>%
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
    select(-.y.) -> geomic_feature_difference_order
geomic_feature_difference_order

addWorksheet(tableS4,sheetName = 'order')
writeData(tableS4, sheet = 'order', x = geomic_feature_difference_order)

## family
df %>%
  filter(Assembly != "GCA_024752465.1") %>%
  filter(Family %in% n_10_family_lst) %>%
  select(Family, Genome_size:tRNA_length_mean, -c(Contig_scaffold_number, N50, BUSCO, Protein_number, GO_number)) %>%
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
    select(-.y.)-> geomic_feature_difference_family
geomic_feature_difference_family


addWorksheet(tableS4,sheetName = 'family',)
writeData(tableS4, sheet = 'family', x = geomic_feature_difference_family, )

## lifestyle
df %>%
  filter(Assembly != "GCA_024752465.1", Lifestyle != "undetermined") %>%
  filter(Lifestyle %in% n_10_lifestyle_lst) %>%
  select(Lifestyle, Genome_size:tRNA_length_mean, -c(Contig_scaffold_number, N50, BUSCO, Protein_number, GO_number)) %>%
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
    select(-.y.) -> geomic_feature_difference_lifestyle
geomic_feature_difference_lifestyle

addWorksheet(tableS4,sheetName = 'lifestyle')
writeData(tableS4, sheet = 'lifestyle', x = geomic_feature_difference_lifestyle)


# dunntest
## subclass
df %>%
  filter(Assembly != "GCA_024752465.1") %>%
  filter(Subclass %in% n_10_subclass_lst) %>%
  select(Subclass, Genome_size:tRNA_length_mean, -c(Contig_scaffold_number, N50, BUSCO, Protein_number, GO_number)) %>%
  pivot_longer(-Subclass, names_to = "feature", values_to = "value") %>%
  group_by(feature) %>%
  dunn_test(value ~ Subclass) %>%
  as.data.frame() %>%
  mutate(different = ifelse(p.adj <= 0.05, "yes", "no")) %>% 
    select(-.y.) -> geomic_feature_pair_difference_subclass
geomic_feature_pair_difference_subclass

df %>%
    filter(Assembly != "GCA_024752465.1") %>%
    filter(Subclass %in% n_10_subclass_lst) %>%
    select(Subclass, Genome_size:tRNA_length_mean, 
           -c(Contig_scaffold_number, N50, BUSCO, Protein_number, GO_number)) %>%
    pivot_longer(-Subclass, names_to = "feature", values_to = "value") %>%
    group_by(Subclass,feature) %>% 
    summarise(average = mean(value)) -> geomic_feature_average_subclass

geomic_feature_average_subclass

geomic_feature_pair_difference_subclass %>% 
    left_join(geomic_feature_average_subclass, by = c('group1' = 'Subclass', 'feature')) %>%
    rename(group1_average=average) %>% 
    left_join(geomic_feature_average_subclass, by = c('group2' = 'Subclass', 'feature')) %>%
    rename(group2_average=average) -> geomic_feature_pair_difference_subclass

geomic_feature_pair_difference_subclass %>% 
    group_by(group1,group2) %>% 
    count(different) %>% arrange(desc(different), desc(n)) -> difference_class2_class
difference_class2_class

difference_class2_class
addWorksheet(tableS4,sheetName = 'class-class')
writeData(tableS4, sheet = 'class-class', x = difference_class2_class)

addWorksheet(tableS4,sheetName = 'pairwise-subclass')
writeData(tableS4, sheet = 'pairwise-subclass', x = geomic_feature_pair_difference_subclass)
geomic_feature_pair_difference_subclass

addWorksheet(tableS4,sheetName = 'pairwise-subclass')
writeData(tableS4, sheet = 'pairwise-subclass', x = geomic_feature_pair_difference_subclass)


## order
df %>%
  filter(Assembly != "GCA_024752465.1") %>%
  filter(Order %in% n_10_order_lst) %>%
  select(Order, Genome_size:tRNA_length_mean, -c(Contig_scaffold_number, N50, BUSCO, Protein_number, GO_number)) %>%
  pivot_longer(-Order, names_to = "feature", values_to = "value") %>%
  group_by(feature) %>%
  dunn_test(value ~ Order) %>%
  as.data.frame() %>%
  mutate(different = ifelse(p.adj <= 0.05, "yes", "no")) %>% 
    select(-.y.) -> genomic_feature_pair_difference_order
n_10_order_lst
genomic_feature_pair_difference_order

df %>%
    filter(Assembly != "GCA_024752465.1") %>%
    filter(Order %in% n_10_order_lst) %>%
    select(Order, Genome_size:tRNA_length_mean, 
           -c(Contig_scaffold_number, N50, BUSCO, Protein_number, GO_number)) %>%
    pivot_longer(-Order, names_to = "feature", values_to = "value") %>%
    group_by(Order,feature) %>% 
    summarise(average = mean(value)) -> genomic_feature_average_order

genomic_feature_average_order

genomic_feature_pair_difference_order %>% 
    left_join(genomic_feature_average_order, by = c('group1' = 'Order', 'feature')) %>%
    rename(group1_average=average) %>% 
    left_join(genomic_feature_average_order, by = c('group2' = 'Order', 'feature')) %>%
    rename(group2_average=average) -> genomic_feature_pair_difference_order
genomic_feature_pair_difference_order

addWorksheet(tableS4,sheetName = 'pairwise-order')
writeData(tableS4, sheet = 'pairwise-order', x = genomic_feature_pair_difference_order)

genomic_feature_pair_difference_order %>% group_by(group1,group2) %>% 
    count(different) %>% arrange(desc(different), desc(n)) -> difference_order2_order
difference_order2_order

addWorksheet(tableS4,sheetName = 'order-order')
writeData(tableS4, sheet = 'order-order', x = difference_order2_order)


## family
df %>%
  filter(Assembly != "GCA_024752465.1") %>%
  filter(Family %in% n_10_family_lst) %>%
  select(Family, Genome_size:tRNA_length_mean, -c(Contig_scaffold_number, N50, BUSCO, Protein_number, GO_number)) %>%
  pivot_longer(-Family, names_to = "feature", values_to = "value") %>%
  group_by(feature) %>%
  dunn_test(value ~ Family) %>%
  as.data.frame() %>%
  mutate(different = ifelse(p.adj <= 0.05, "yes", "no")) %>% 
    select(-.y.) -> genomic_feature_pair_difference_family
n_10_family_lst
genomic_feature_pair_difference_family

df %>%
    filter(Assembly != "GCA_024752465.1") %>%
    filter(Family %in% n_10_family_lst) %>%
    select(Family, Genome_size:tRNA_length_mean, 
           -c(Contig_scaffold_number, N50, BUSCO, Protein_number, GO_number)) %>%
    pivot_longer(-Family, names_to = "feature", values_to = "value") %>%
    group_by(Family,feature) %>% 
    summarise(average = mean(value)) -> genomic_feature_average_family

genomic_feature_average_family

genomic_feature_pair_difference_family %>% 
    left_join(genomic_feature_average_family, by = c('group1' = 'Family', 'feature')) %>%
    rename(group1_average=average) %>% 
    left_join(genomic_feature_average_family, by = c('group2' = 'Family', 'feature')) %>%
    rename(group2_average=average) -> genomic_feature_pair_difference_family
genomic_feature_pair_difference_family


addWorksheet(tableS4,sheetName = 'pairwise-family')
writeData(tableS4, sheet = 'pairwise-family', x = genomic_feature_pair_difference_family)

genomic_feature_pair_difference_family %>% group_by(group1,group2) %>% 
    count(different) %>% arrange(desc(different), desc(n)) -> difference_family2_family
difference_family2_family

addWorksheet(tableS4,sheetName = 'family-family')
writeData(tableS4, sheet = 'family-family', x = difference_family2_family)


## lifestyle

df %>%
  filter(Assembly != "GCA_024752465.1", Lifestyle != "undetermined") %>%
  filter(Lifestyle %in% n_10_lifestyle_lst) %>%
  select(Lifestyle, Genome_size:tRNA_length_mean, -c(Contig_scaffold_number, N50, BUSCO, Protein_number, GO_number)) %>%
  pivot_longer(-Lifestyle, names_to = "feature", values_to = "value") %>%
  group_by(feature) %>%
  dunn_test(value ~ Lifestyle) %>%
  as.data.frame() %>%
  mutate(difference = ifelse(p.adj <= 0.05, "yes", "no")) %>% 
    select(-.y.) -> genomic_feature_pair_difference_lifestyle
genomic_feature_pair_difference_lifestyle

df %>%
    filter(Assembly != "GCA_024752465.1") %>%
    filter(Lifestyle %in% n_10_lifestyle_lst) %>%
    select(Lifestyle, Genome_size:tRNA_length_mean, 
           -c(Contig_scaffold_number, N50, BUSCO, Protein_number, GO_number)) %>%
    pivot_longer(-Lifestyle, names_to = "feature", values_to = "value") %>%
    group_by(Lifestyle,feature) %>% 
    summarise(average = mean(value)) -> genomic_feature_average_lifestyle

genomic_feature_average_lifestyle

genomic_feature_pair_difference_lifestyle %>% 
    left_join(genomic_feature_average_lifestyle, by = c('group1' = 'Lifestyle', 'feature')) %>%
    rename(group1_average=average) %>% 
    left_join(genomic_feature_average_lifestyle, by = c('group2' = 'Lifestyle', 'feature')) %>%
    rename(group2_average=average) -> genomic_feature_pair_difference_lifestyle
genomic_feature_pair_difference_lifestyle

addWorksheet(tableS4,sheetName = 'pairwise-lifestyle')
writeData(tableS4, sheet = 'pairwise-lifestyle', x = genomic_feature_pair_difference_lifestyle)

genomic_feature_pair_difference_lifestyle  %>% group_by(group1,group2) %>% 
    count(difference) %>% arrange(desc(difference), desc(n)) -> difference_lifestyle2_lifestyle
difference_lifestyle2_lifestyle
addWorksheet(tableS4,sheetName = 'lifestyle-lifestyle')
writeData(tableS4, sheet = 'lifestyle-lifestyle', x = difference_lifestyle2_lifestyle)


saveWorkbook(tableS4, "TableS4_genomic_feature_difference.xlsx", overwrite = TRUE)

# plot
## subclass
feature_lst <- rev(c(
    "Genome_size",  "Genome_size_without_TE", 'TE_size',
    "GC_content", "GC_content_without_TE", "GC_content_TE",
    "tRNA_length_max", "tRNA_length_mean", "tRNA_length_min", "tRNA_number",
    "Intergenic_length_max", "Intergenic_length_mean", "Intergenic_length_min",
    "Intron_length_max", "Intron_length_mean", "Intron_length_min", "Intron_number",
    "Exon_length_max", "Exon_length_mean", "Exon_length_min","Exon_number",
    "Gene_length_max", "Gene_length_mean", "Gene_length_min", "Gene_number"
))

geomic_feature_pair_difference_subclass




geomic_feature_pair_difference_subclass %>%
    group_by(feature, different) %>%
    count() %>%
    as.data.frame() %>% 
    mutate(level='subclass') -> geomic_feature_pair_difference_subclass_count_df

geomic_feature_pair_difference_subclass_count_df

geomic_feature_pair_difference_subclass %>%
  group_by(feature, different) %>%
  count() %>%
  as.data.frame() %>%
  ggplot(.) +
  geom_bar(aes(y = factor(feature, levels = rev(cluster_order_lst)), 
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
geomic_feature_pair_difference_order %>%
    group_by(feature, different) %>%
    count() %>%
    as.data.frame() %>%
    mutate(level='order') -> geomic_feature_pair_difference_order_count_df
geomic_feature_pair_difference_order_count_df

geomic_feature_pair_difference_order %>%
  group_by(feature, different) %>%
  count() %>%
  as.data.frame() %>%
  ggplot(.) +
  geom_bar(aes(y = factor(feature, levels = rev(cluster_order_lst)), 
               x = n, fill = different), 
           position = "stack", 
           stat = "identity") +
  theme_minimal() + 
    ylab(NULL) + 
    theme(axis.text.y = element_blank()) + 
    xlab(NULL) + ggtitle('Order') -> pair_difference_order_p
pair_difference_order_p

## family
geomic_feature_pair_difference_family %>%
    group_by(feature, different) %>%
    count() %>%
    as.data.frame() %>% 
    mutate(level='family') -> geomic_feature_pair_difference_family_count_df
geomic_feature_pair_difference_family_count_df

geomic_feature_pair_difference_family %>%
  group_by(feature, different) %>%
  count() %>%
  as.data.frame() %>%
  ggplot(.) +
  geom_bar(aes(y = factor(feature, levels = rev(cluster_order_lst)), 
               x = n, 
               fill = different), 
           position = "stack", 
           stat = "identity") +
  theme_minimal() + ylab(NULL) +
    theme(axis.text.y = element_blank()) + 
    xlab(NULL) + ggtitle('Family') -> pair_difference_family_p
pair_difference_family_p


## lifestyle
geomic_feature_pair_difference_lifestyle %>%
    group_by(feature, difference) %>%
    count() %>%
    as.data.frame() %>%
    mutate(level='lifestyle') %>% 
    rename(different=difference)-> geomic_feature_pair_difference_lifestyle_count_df
geomic_feature_pair_difference_lifestyle_count_df

geomic_feature_pair_difference_lifestyle %>%
  group_by(feature, difference) %>%
  count() %>%
  as.data.frame() %>%
  ggplot(.) +
  geom_bar(aes(y = factor(feature, levels = rev(cluster_order_lst)), 
               x = n, 
               fill = difference), 
           position = "stack", 
           stat = "identity") +
  theme_minimal() + ylab(NULL) + 
    theme(axis.text.y = element_blank()) + 
    xlab(NULL) + ggtitle('Lifestyle')-> pair_difference_lifestyle_p
pair_difference_lifestyle_p

## merge

rbind(geomic_feature_pair_difference_subclass_count_df, 
      geomic_feature_pair_difference_order_count_df,
      geomic_feature_pair_difference_family_count_df,
      geomic_feature_pair_difference_lifestyle_count_df) -> merge_df

merge_df %>% pivot_wider(names_from = c('level', 'different'), values_from = 'n') %>%
    mutate_at(c(2:9), ~replace_na(.,0)) %>% column_to_rownames('feature') -> significant_count_df

addWorksheet(tableS4,sheetName = 'clustering-matrix')
writeData(tableS4, sheet = 'clustering-matrix', x = significant_count_df)


saveWorkbook(tableS4, "TableS4_genomic_feature_difference.xlsx", overwrite = TRUE)


merge_df %>% pivot_wider(names_from = c('level', 'different'), values_from = 'n') %>%
    mutate_at(c(2:9), ~replace_na(.,0)) %>% column_to_rownames('feature') %>%
    scale() %>% dist(method = 'euclidean')  %>% 
        hclust() %>% as.phylo() %>% plot() -> cluster_p
cluster_p

help(hclust)
help(dist)
help(plot)

cluster_order_lst <- c("tRNA_length_mean", "Exon_length_max", "Intron_number",
                       "Genome_size", "GC_content_without_TE",  "GC_content", "Exon_length_mean",
                       "Gene_number", "Exon_number", "Intron_length_mean", 
                       "Gene_length_mean", "tRNA_length_max", "Gene_length_max",
                       "GC_content_TE", "tRNA_number", "Genome_size_without_TE",
                       "tRNA_length_min", 'TE_size', "Gene_length_min",  "Exon_length_min",
                       "Intron_length_min", "Intergenic_length_max", "Intron_length_max",
                       "Intergenic_length_min", "Intergenic_length_mean")


ggarrange(pair_difference_subclass_p, pair_difference_order_p, 
          pair_difference_family_p, pair_difference_lifestyle_p, 
          nrow = 1, common.legend = TRUE, widths = c(1,1,1,1)) 
