setwd("G:\\Experiments\\Paper\\FD_lifestyle_changes\\11_stat\\machine_learning")
dir()

library(ggplot2)
library(ggsci)
library(tidyverse)
library(pheatmap)
library(ggpubr)
library(rstatix)
library(scales)
library(openxlsx)

# read data table
df1 <- read.csv("classifiers_accuracy_basic_genomic_features.csv")
df1['dataset'] <- 'genomic_feature'
df1
df2 <- read.csv("classifiers_accuracy_function_groups.csv")
df2['dataset'] <- 'functional_protein'
df2
df3 <- read.csv("classifiers_accuracy_together.csv")
df3['dataset'] <- 'Joined'
df3

tableS7 <- createWorkbook()
addWorksheet(tableS7, sheetName = 'accuracy_genomic_feature')
writeData(tableS7, sheet = 'accuracy_genomic_feature', x= df1)

addWorksheet(tableS7, sheetName = 'accuracy_functional_protein')
writeData(tableS7, sheet = 'accuracy_functional_protein', x= df2)

addWorksheet(tableS7, sheetName = 'accuracy_joined')
writeData(tableS7, sheet = 'accuracy_joined', x= df3)

merged_df <- rbind.data.frame(df1,df2,df3)
merged_df %>% group_by(dataset, Classifier) %>% 
    shapiro_test(Accuracy) %>%
    as.data.frame() %>%
    mutate(significance = case_when(
        p <= 0.0001 ~ "****",
        p <= 0.001 ~ "***",
        p <= 0.01 ~ "**",
        p <= 0.05 ~ "*",
        p > 0.05 ~ "ns"
    ),normality = ifelse(p > 0.05, "yes", "no")) %>%
    rename(classifier=Classifier) -> normality_test_df

normality_test_df

addWorksheet(tableS7, sheetName = 'normality-test')
writeData(tableS7, sheet = 'normality-test', x= normality_test_df)

## dunntest
merged_df %>% group_by(dataset, Classifier) %>% 
    summarise(average=mean(Accuracy)) -> average_accuracy_df
average_accuracy_df
addWorksheet(tableS7, sheetName = 'average-accuracy')
writeData(tableS7, sheet = 'average-accuracy', x= average_accuracy_df)


merged_df %>% group_by(dataset) %>% 
    dunn_test(Accuracy ~ Classifier) %>% as.data.frame() %>% 
    mutate(different = ifelse(p.adj > 0.05, 'no', 'yes')) %>% 
    left_join(average_accuracy_df, by = c('dataset', 'group1' = 'Classifier')) %>%
    rename(group1_mean=average) %>%
    left_join(average_accuracy_df, by = c('dataset', 'group2' = 'Classifier')) %>%
    rename(group2_mean=average) %>% select(-.y.) -> dunntest_method_df
dunntest_method_df
addWorksheet(tableS7, sheetName = 'method-method')
writeData(tableS7, sheet = 'method-method', x= dunntest_method_df)

merged_df %>% group_by(Classifier) %>% 
    dunn_test(Accuracy ~ dataset) %>% as.data.frame() %>% 
    mutate(different = ifelse(p.adj > 0.05, 'no', 'yes')) %>% 
    left_join(average_accuracy_df, by = c('Classifier', 'group1' = 'dataset')) %>%
    rename(group1_mean=average) %>%
    left_join(average_accuracy_df, by = c('Classifier', 'group2' = 'dataset')) %>%
    rename(group2_mean=average) %>% select(-.y.) -> dunntest_dataset_df
dunntest_dataset_df
addWorksheet(tableS7, sheetName = 'dataset-dataset')
writeData(tableS7, sheet = 'dataset-dataset', x= dunntest_dataset_df)
saveWorkbook(tableS7, file = 'TableS7_machine_learning_accuracies.xlsx', overwrite = T)


## supplementary figure

merged_df %>% select(-X) %>%
    ggplot(aes(x=dataset, y=Accuracy)) + 
    geom_boxplot(aes(color=dataset)) + 
    facet_wrap(~Classifier) + theme_bw() +
    theme(axis.text.x = element_text(angle=45, vjust = 0.5)) + 
    scale_color_lancet()
    
    


# plot
ggboxplot(df1, x = "Classifier", y = "Accuracy",
          color = "Classifier", add = "jitter") + 
    scale_color_lancet() + 
    theme_minimal() + 
    geom_hline(yintercept = mean(df1$Accuracy), linetype = 2) + 
    geom_label(x=1, y = mean(df1$Accuracy), label = round(mean(df1$Accuracy),4)) + 
    ggtitle('Based on basic genomic features') + 
    theme(legend.position = 'none') -> p1
p1

ggboxplot(df2, x = "Classifier", y = "Accuracy",
          color = "Classifier", add = "jitter") + 
    scale_color_lancet() + 
    theme_minimal() + 
    geom_hline(yintercept = mean(df2$Accuracy), linetype = 2) + 
    geom_label(x=1, y = mean(df2$Accuracy), label = round(mean(df2$Accuracy),4)) + 
    ggtitle('Based on functional protein groups') +
    theme(legend.position = 'none') -> p2
p2


ggboxplot(df3, x = "Classifier", y = "Accuracy",
          color = "Classifier", add = "jitter") + 
    scale_color_lancet() + 
    theme_minimal() + 
    geom_hline(yintercept = mean(df3$Accuracy), linetype = 2) + 
    geom_label(x=1, y = mean(df3$Accuracy), label = round(mean(df3$Accuracy),4)) + 
    ggtitle('Joined dataset')  + 
    theme(legend.position = 'none')-> p3
p3

# for legend
ggboxplot(df3, x = "Classifier", y = "Accuracy",
          color = "Classifier", add = "jitter") + 
    scale_color_lancet() + 
    theme_minimal() + 
    geom_hline(yintercept = mean(df3$Accuracy), linetype = 2) + 
    geom_label(x=1, y = 0.7287, label = round(mean(df3$Accuracy),4)) + 
    ggtitle('Joined dataset') -> p3_legend
    


ggarrange(p1,p2,p3, ncol = 1) -> pm1
pm1
## confusion matrix
dc1 <- read.csv("confusion_matrix_basic_genomic_features.csv", header = T, row.names = 1)
dc2 <- read.csv("confusion_matrix_functional_protein_groups.csv", header = T, row.names = 1)
dc3 <- read.csv("confusion_matrix_together.csv", header = T, row.names = 1)

dc1 %>%
  mutate(newSum = select_if(., is.numeric) %>%
    reduce(`+`)) %>%
  mutate_if(is.numeric, list(~ . / newSum * 100)) %>%
  select(-newSum) %>%
  pheatmap(., cluster_rows = F, cluster_cols = F, 
           border = "white", 
           angle_col = 45, 
           display_numbers = T, 
           number_format = "%.2f", legend = F)  -> pc1
pc1

help(pheatmap)
dc2 %>%
    mutate(newSum = select_if(., is.numeric) %>%
               reduce(`+`)) %>%
    mutate_if(is.numeric, list(~ . / newSum * 100)) %>%
    select(-newSum) %>%
    pheatmap(., cluster_rows = F, cluster_cols = F, 
             border = "white", 
             angle_col = 45, 
             display_numbers = T, 
             number_format = "%.2f", 
             show_rownames = T, legend = F) -> pc2
pc2

dc3 %>%
    mutate(newSum = select_if(., is.numeric) %>%
               reduce(`+`)) %>%
    mutate_if(is.numeric, list(~ . / newSum * 100)) %>%
    select(-newSum) %>%
    pheatmap(., cluster_rows = F, cluster_cols = F, 
             border = "white", 
             angle_col = 45, 
             display_numbers = T, 
             number_format = "%.2f", 
             show_rownames = T,
             legend = F) -> pc3
pc3

# for legend
dc3 %>%
    mutate(newSum = select_if(., is.numeric) %>%
               reduce(`+`)) %>%
    mutate_if(is.numeric, list(~ . / newSum * 100)) %>%
    select(-newSum) %>%
    pheatmap(., cluster_rows = F, cluster_cols = F, 
             border = "white", 
             angle_col = 45, 
             display_numbers = T, 
             number_format = "%.2f", 
             show_rownames = T,
             legend = T) -> pc3_legend

ggarrange(pc1$gtable, pc2$gtable, pc3$gtable,ncol = 1) -> pm2

ggarrange(p1, pc1$gtable, p2, pc2$gtable,p3, pc3$gtable, 
          ncol = 2, nrow = 3, 
          widths = c(1,1,1,1,1,1), 
          labels = c('a','b','c','d','e','f'))

ggarrange(p3_legend, pc3_legend$gtable, ncol = 2)


##predict undetermined
predict <- read.xlsx('TableS8_predicted_lifestyle.xlsx')
predict %>% group_by(predicted_lifestyle, Subclass, Order, Family) %>% 
    count() %>% as.data.frame()


