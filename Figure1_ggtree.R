# set working directory
setwd('G:\\Experiments\\Paper\\FD_lifestyle_changes\\11_stat')

# list content in the directory
dir()

# load required packages
library(treedataverse)
library(tidyverse)
library(readxl)
library(ggstar)
library(ggnewscale)

# read in the ML tree
tre <- read.iqtree('G:\\Experiments\\Paper\\FD_lifestyle_changes\\12_ml_tree\\640_1124_iqtree.treefile')

# read in the meadadata table
df <- read_excel('G:\\Experiments\\Paper\\FD_lifestyle_changes\\ncbi_genome_table_20220812\\20221107.xlsx')

# convert the sequencing technologies into factor variables
df$WGS_class <- factor(df$WGS_class)

# to label clades
highlight_node_lst <- c(647,1010,1051,1116,1157,1167,1205,1226,1243,1229)
node_label <- c('Hypocreales','Microascales', 'Glomerellales','Ophiostomatales','Magnaporthales',
                'Diaporthales', 'Sordariales', 'Coniochaetales','Xylariales','Amphisphaeriales')
nodedf <- data.frame(node=highlight_node_lst,label=node_label)

# tree + calde labels
p <- ggtree(tre, layout = 'fan',size=0.15) +
    geom_cladelab(data=nodedf, 
                  mapping=aes(node=node, 
                              label=label), fontsize=2, offset.text=0.1)
p


# terminal shapes represent the lifestyles
p1 <- p + geom_fruit(data = df, geom = geom_star,
               mapping = aes(y=Assembly, fill=Lifestyle, starshape=WGS_class),
               position='identity',
               starstroke=0.2) + theme(legend.position = 'bottom')

df %>% select(Assembly,Gender) %>% 
    mutate(sexual = case_when(Gender == 'sexual' ~ 1,
                              Gender == 'sexual & asexual' ~ 1,
                              TRUE ~ 0)) %>% 
    mutate(asexual = case_when(Gender == 'asexual' ~ 1,
                              Gender == 'sexual & asexual' ~ 1,
                              TRUE ~ 0)) %>% select(Assembly,sexual,asexual) %>% 
    gather(key='Gender', value = 'Values', -Assembly) -> gender_df
typeof(gender_df$Values)

### we already removethe morph information frome the figure
p2 <- p1 + geom_fruit(data = gender_df,geom = geom_tile,
                mapping = aes(y = Assembly, x = Gender, alpha=Values),
                pwidth=0.02, color = "grey50", offset = 0.04,size = 0.02, fill='steelblue')  + 
    scale_alpha_continuous(name='Sexual stage',range = c(0,1))

# gc content line
df %>% select(Assembly, GC_content,GC_content_without_TE) %>% 
    gather(key='withTE', value = 'content',-Assembly) %>%
    mutate(content = replace(content, is.na(content), '30'))  -> df_gc

range(df_gc$content)


p3 <- p2 + geom_fruit(data = df_gc, geom_line,
                mapping = aes(y=Assembly, x=content,color=withTE), position=position_identityx()) 

p3
# genome size
df %>% select(Assembly, Genome_size, TE_size) %>% 
    mutate(genome_size_no_te=Genome_size-TE_size) %>% 
    select(Assembly, genome_size_no_te,TE_size) %>%
    gather(key = 'type', value = 'size', -Assembly)-> df_genome

range(df_genome$size)

df_genome$type <- factor(df_genome$type, levels = c('TE_size','genome_size_no_te' ))

df_genome$type

p3 + geom_fruit(data = df_genome, geom = geom_bar,
                mapping = aes(y=Assembly, x=size, color=type), pwidth=0.1,fill='gray', 
                orientation="y", 
                stat="identity") + 
    geom_treescale(fontsize=2, linesize=0.5,x=0,y=0)

