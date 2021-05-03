#################################################################################
#################################################################################
################## Monocle3 (Tweets from four distinct users) ###################
#################################################################################
#################################################################################

setwd("C:/Users/jueju/Desktop/LDA and Co-clustering algorithms with data")
getwd()

library(dplyr)
library(monocle3)
library(cluster)
library(ggplot2)

####### expression ########
sparse_matrix = read.csv("doc_word_matrix_stemmingf.csv") %>% dplyr::select(-X)
names = names(sparse_matrix)
sparse_matrix = t(as.matrix(sparse_matrix))
word_expression = sparse_matrix

dim(word_expression)
word_expression[1:5, 1:5]

####### Doc_metadata #######
doc_metadata = read.csv("df_doc_topic_new_stemming.csv")
doc_metadata = doc_metadata  %>% dplyr::select(Keywords, Doc_Text, screen_name, user_id_new)
dim(doc_metadata)
head(doc_metadata)

####### Word_metadata #######
word_metadata = data.frame(
  word = names,
  gene_short_name = names
)

dim(word_metadata)
head(word_metadata)

######## Change row_col names#######
colnames(word_expression) = rownames(doc_metadata)
rownames(word_metadata) = rownames(word_expression)

######## Create Monocle3 Object ##########
cds = new_cell_data_set(expression_data = word_expression,
                         cell_metadata = doc_metadata,
                         gene_metadata = word_metadata)

######## Pre-process ##########
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds)

#####################################
############ scree plot #############
#####################################

Monocle3_PCA = cds@int_colData@listData$reducedDims@listData$PCA
nnn = as.data.frame(Monocle3_PCA)
vars_transformed <- apply(nnn, 2, var)
Monocle3_PCA = vars_transformed/sum(vars_transformed)
Monocle3_PCA = as.numeric(Monocle3_PCA)

PCA <- data.frame(
  PCA_components = 1:100,
  Monocle3 = Monocle3_PCA
)

PCA %>% ggplot(aes(x = PCA_components, y = Monocle3)) + geom_point() +
  ylab("Proportion of variance explained by each component") +
  ggtitle("Monocle3")

#################################################
############# Silhouette Analysis ###############
#################################################

############ Distance matrix ############## 

# (1) Choice 1: Ferg's distance matrix
distance_matrix <- read.csv("Robyn_tweets_distance.csv", header = FALSE)

# (2) Choice 2: Euclidean distance matrix
reduced_dim_res <- reducedDims(cds)[["UMAP"]]
distance_matrix = dist(reduced_dim_res, diag = T, upper = T)

############ iterations ############
k_range  =c(5,10,20,30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,240,260,280)
avg_widths = vector()


for (i in 1:length(k_range)){
  num = as.integer(k_range[i])
  cds = cluster_cells(cds,  k = num)
  x = cds@clusters@listData$UMAP$cluster_result$optim_res$membership
  si = silhouette(x, distance_matrix)
  ssi <- summary(si)
  avg_widths[i] = ssi$avg.width
}

library(ggplot2)
summary_s = data.frame(k_range, avg_widths)
str(summary_s)
ggplot(data=summary_s, aes(x=k_range, y=avg_widths)) +
  geom_line()+
  geom_point() + scale_x_continuous(name="k", breaks = k_range) +
  scale_y_continuous(name="Average Silhouette Width") + ggtitle("Silhouette for Monocle3")

########################################################
############# Useful Visualization Tools ###############
########################################################

k = k_range[which.max(avg_widths)]
cds = cluster_cells(cds,  k)

### UMAP ###
# Clustering results
plot_cells(cds, group_label_size = 5)

# True labels
plot_cells(cds, color_cells_by="screen_name", group_label_size = 5, label_groups_by_cluster = FALSE)

### Feature plot ###
plot_cells(cds, genes=c("storm", "health", "goblu", "vegan"))


