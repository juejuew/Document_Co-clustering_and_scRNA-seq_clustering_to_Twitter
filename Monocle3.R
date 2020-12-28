###############################################
###############################################
################## Monocle3 ###################
############################################### 
############################################### 

setwd("C:/Users/jueju/Desktop/LDA and Co-clustering algorithms with data")
getwd()

library(dplyr)
library(monocle3)
library(cluster)

####### expression ########
sparse_matrix <- read.csv("test_new_stem.csv", header = FALSE)
sparse_matrix = as.matrix(sparse_matrix)
sparse_matrix = t(sparse_matrix)
word_expression = sparse_matrix
dim(word_expression)

####### Doc_metadata
doc_metadata <- read.csv("df_doc_topic_new_stemming.csv")
dim(doc_metadata)
str(doc_metadata)
doc_metadata  = doc_metadata  %>% select(Keywords, Doc_Text, screen_name, user_id_new)
dim(doc_metadata)

####### Word_metadata
word_metadata <- read.csv("word_metadata_new.csv")
word_metadata = word_metadata %>% select(word) %>% mutate(gene_short_name = word)
str(word_metadata)
dim(word_metadata)

######## Change row_col names#######
rownames(word_expression) <- rownames(word_metadata)
colnames(word_expression) <- rownames(doc_metadata)

######## Pre-process ##########
cds <- new_cell_data_set(expression_data = word_expression,
                         cell_metadata = doc_metadata,
                         gene_metadata = word_metadata)

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
distance_matrix <- read.csv("Robyn_tweets_distance.csv", header = FALSE)

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

max(avg_widths)
length(k_range)
length(avg_widths)

library(ggplot2)
omg = data.frame(k_range, avg_widths)
ggplot(data=omg, aes(x=k_range, y=avg_widths)) +
  geom_line()+
  geom_point() + scale_x_continuous(name="k", breaks = k_range) +
  scale_y_continuous(name="Average Silhouette Width") + ggtitle("Silhouette for Monocle3")

########################################################
############# Useful Visualization Tools ###############
########################################################

cds = cluster_cells(cds,  k = 65)

# UMAP
plot_cells(cds)
plot_cells(cds, color_cells_by="screen_name")

# Marker words by each cluster
marker_test_res <- top_markers(cds, 
                               reference_cells=1000, cores=8)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(5, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    #group_cells_by="partition",
                    ordering_type="cluster_row_col",
                    max.size=5)
