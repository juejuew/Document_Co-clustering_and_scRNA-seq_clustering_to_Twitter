###############################################################
###############################################################
############# Proposed Workflow ("jobs" tweets) ###############
###############################################################
###############################################################


###############################
######## Combine data #########
###############################

library(dplyr)
library(tidyverse)
library(data.table)
library(lubridate)
library(ggplot2)
library(stringr)
library(monocle3)
library(cluster)

# Sampling
setwd("C:/Users/jueju/Desktop/jobs tweets/ByDay")
getwd()

df = list.files(pattern = ".rds") %>%
  map(readRDS) %>% 
  data.table::rbindlist()

df_sample_1 = df %>% dplyr::filter(time < "2009-10-31 23:59:59" & time > "2009-7-31 23:59:59") %>%
  mutate(date = as.Date(ymd_hms(time)))

########################
####### Cleaning #######
########################

####### Remove retweets ########

no.rts = df_sample_1[grep("^RT ", df_sample_1$text, invert=TRUE),]

####### Remove duplicates ########

df_after_cleaning = no.rts[!duplicated(no.rts$text),]

# Plot
df_after_cleaning %>% ggplot +
  geom_bar(aes(x = date)) +
  scale_y_continuous(limits=c(0, 15000))

###########################
######## Sampling #########
###########################

set.seed(88)
df_sample = df_after_cleaning %>% group_by(date) %>% sample_n(300) %>% as.data.frame(.) #%>% select(-n)

# Plot
df_sample %>% ggplot +
  geom_bar(aes(x = date)) 

#write.csv(Target,"f_jobs_tweets_sampled_three_month.csv")


##############################################################################################################################################################################################
#Preprocessing can be found on https://github.com/juejuew/Document_Co-clustering_and_scRNA-seq_clustering_to_Twitter/blob/main/proposed_workflow/Proposed%20workflow%20-%20preprocessing.ipynb
##############################################################################################################################################################################################

setwd("C:/Users/jueju/Desktop/jobs tweets")
getwd()

####### expression ########
sparse_matrix = read.csv("f_jobs_doc_word_matrix_stemmingf.csv")
sparse_matrix = sparse_matrix %>% select(-X)
names = names(sparse_matrix)
sparse_matrix = t(as.matrix(sparse_matrix))
word_expression = sparse_matrix
dim(word_expression)

####### Word_metadata
word_metadata = data.frame(
  word = names,
  gene_short_name = names
)

dim(word_metadata)

####### Doc_metadata
doc_metadata = read.csv("f_jobs_stemming_meta_doc.csv") %>% mutate(label = as.character(label))
dim(doc_metadata)

####### Create the Monocle3 object
colnames(word_expression) = rownames(doc_metadata)
rownames(word_metadata) = rownames(word_expression)

cds = new_cell_data_set(expression_data = word_expression,
                         cell_metadata = doc_metadata,
                         gene_metadata = word_metadata)


cds = preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)
cds = reduce_dimension(cds)

##########################
###### UMAP plot #########
##########################

plot_cells(cds)
plot_cells(cds, color_cells_by="label", group_label_size = 5, label_groups_by_cluster = FALSE)
 
###### Study the content of the small communities
selected_tweets_idx = choose_cells(cds, return_list = TRUE)
selected_tweets_info = doc_metadata %>% mutate(row_num =row_number()) %>% filter(row_num %in% as.numeric(selected_tweets_idx))

##########################
######### DBSCAN #########
##########################

###### Apply DBSCAN to seperate and store the main cluster and the small communities
reduced_dim_res = reducedDims(cds)[["UMAP"]]

as.data.frame(reduced_dim_res) %>% ggplot + geom_point(aes(x = V1, y = V2), size = 0.1)

#write.csv(as.data.frame(reduced_dim_res), "f_UMAP_idx.csv")

#####################################################################################################################################################################
# Find the implementation of DBSCAN on https://github.com/juejuew/Document_Co-clustering_and_scRNA-seq_clustering_to_Twitter/blob/main/proposed_workflow/DBSCAN.ipynb
#####################################################################################################################################################################

DBSCAN = read.csv("f_label_DBSCAN.csv")
str(DBSCAN)

useful1 = as.data.frame(reduced_dim_res) %>% mutate(X = DBSCAN$X,label_DBSCAN = as.character(DBSCAN$label_DBSCAN))

useful1 %>% 
  ggplot(aes(x = V1, y = V2, color = label_DBSCAN)) + 
  geom_point(size = 0.1, show.legend = FALSE) + 
  geom_text(
    label=useful1$label_DBSCAN, 
    nudge_x = 0.25, nudge_y = 0.25, 
    check_overlap = T
  )  + theme(legend.position="none")


useful2 = useful1 %>% dplyr::count(label_DBSCAN) 
useful2 %>% ggplot + geom_bar(aes(x = reorder(label_DBSCAN, -n), y = n),stat="identity") + coord_flip()

##################################
########## Main Cluster ##########
##################################
select_cluster = sort(useful2$n, decreasing = TRUE)[1:2]
select_label = useful2 %>% filter(n %in% select_cluster) %>% .$label_DBSCAN

main_Finally = useful1 %>% filter(label_DBSCAN %in% select_label)
main_Finally %>% mutate(label_DBSCAN = as.character(label_DBSCAN)) %>%
  ggplot + geom_point(aes(x = V1, y = V2, color = label_DBSCAN), size = 0.1, show.legend = FALSE)

main_Finally = main_Finally %>% mutate(X = X + 1)
str(main_Finally)
main = cds[,main_Finally$X]

dim(main)

#write.csv(doc_metadata[main_Finally$X,], "fsub_jobs_main_meta_doc_stemming.csv")
#write.csv(t(sparse_matrix)[main_Finally$X,], "fsub_jobs_main_doc_word_matrix_stemming.csv")

#######################################
########## Small communities ##########
#######################################

select_cluster_small = sort(useful2$n, decreasing = TRUE)[3:162]
select_label_small = useful2 %>% filter(n %in% select_cluster_small) %>% .$label_DBSCAN

small_Finally = useful1 %>% filter(label_DBSCAN %in% select_label_small)
small_Finally %>% mutate(label_DBSCAN = as.character(label_DBSCAN)) %>%
  ggplot + geom_point(aes(x = V1, y = V2, color = label_DBSCAN), size = 0.1, show.legend = FALSE)

small_Finally = small_Finally %>% mutate(X = X + 1)
str(small_Finally)
small = cds[,small_Finally$X]

dim(small)

######################################################
######################################################
############ Focusing on the main cluster ############
######################################################
######################################################

##########################################
########### Silhouette Analysis ##########
##########################################

########### distance matrix ############

reduced_dim_res = reducedDims(main)[["UMAP"]]
distance_matrix = dist(reduced_dim_res, diag = T, upper = T)

##########################################

#######################################################################################
# Find the value for the tuning parameter (i.e., the number of nearest neighbors) such that the corresponding number of clusters is 2
set.seed(44)
main = cluster_cells(main,  k = 150)
x = main@clusters@listData$UMAP$cluster_result$optim_res$membership
as.integer(length(unique(x)))

########################################################################################

Monocle3_results = data.frame(1:dim(main)[2])
num_cluster = vector()
avg_widths = vector()

k_range  =c(5,10,20,30,40,50,60,70,80,90,100,110,120,140, 150)


for (i in 1:length(k_range)){
  set.seed(44)
  num = as.integer(k_range[i])
  main = cluster_cells(main,  k = num)
  x = main@clusters@listData$UMAP$cluster_result$optim_res$membership
  
  #Store results
  num_cluster[i] = as.integer(length(unique(x))) # count the number of clusters
  nn = data.frame(as.numeric(x))
  Monocle3_results = cbind(Monocle3_results, nn)
  
  # Silhouette Analysis
  si = silhouette(x, distance_matrix)
  ssi = summary(si)
  avg_widths[i] = ssi$avg.width
}

summary_s  = data.frame(k_range, avg_widths, num_cluster)
names(summary_s) = c("k_range", "avg_widths", "num_cluster")

ggplot(data=summary_s , aes(x=k_range, y=avg_widths)) +
  geom_line()+
  geom_point() + scale_x_continuous(name="k", breaks = k_range) +
  scale_y_continuous(name="Average Silhouette Width") + ggtitle("Silhouette for Monocle3")


#write.csv(summary_s, "monocle3_jobs_silhouette.csv")
#write.csv(Monocle3_results, "monocle3_jobs_results.csv")

#################################################################
##################### Monocle3 Clustering #######################
#################################################################

k = k_range[which.max(avg_widths)]
main = cluster_cells(main,  k, reduction_method = "UMAP")

### UMAP ###
# Clustering results
plot_cells(main, group_label_size = 5)

# True labels
plot_cells(main, color_cells_by="label", group_label_size = 5, label_groups_by_cluster = FALSE)

##################################################################
##################### Feature word selection #####################
##################################################################

classification_results = read.csv("monocle3_jobs_results.csv")[,-1]
names(classification_results) = c("idx", k_range)
df_compare = read.csv("fsub_jobs_main_meta_doc_stemming.csv")
dim(df_compare)

# Contingency table
classification_results_selected = as.numeric(classification_results[,which.max(avg_widths)])
df_compare = df_compare %>% mutate(Monocle3_results = classification_results_selected)
organize = df_compare  %>% dplyr::count(category, Monocle3_results)  
contingency_table = spread(organize, category, n, fill = 0, convert = FALSE)
contingency_table

# Summarize word frequecncy in each cluster
new_sparse_matrix = t(sparse_matrix)
summary_table = data.frame(1:ncol(new_sparse_matrix))

for(i in 1:length(unique(df_compare$Monocle3_results))){
  
  cluster_i = new_sparse_matrix[as.numeric(df_compare %>% filter(Monocle3_results == i) %>% .[,1]),]
  
  words = colSums(cluster_i)
  sort_by_freq = sort(words, decreasing = TRUE)
  
  summary_table = cbind(summary_table, data.frame(names(sort_by_freq)))
}

summary_table = summary_table[-1,-1]
names(summary_table) = as.character(1:length(unique(df_compare$Monocle3_results)))

summary_table[1:26,]

# Select feature words
marker = unique(as.vector(as.matrix(summary_table)[1:10,]))

Final_marker = vector()

for(i in 1:length(marker)){
  # exclude the words with <= 3 letters
  if(nchar(marker[i])>3){
    Final_marker = c(Final_marker, marker[i])
  }
}

length(Final_marker)
Final_marker

#write.csv(data.frame(Final_marker), "jobs_main_final_marker.csv")

#####################################################################
##################### Identification of Topics ######################
#####################################################################

##################################################################
# Source code from: https://github.com/cole-trapnell-lab/monocle3
##################################################################

# Marker words by each cluster
Markers = read.csv("jobs_main_final_marker.csv")
markers = Markers$Final_marker

# Define function
plot_words_by_group = function(
  cds,
  markers,
  opt_k,
  group_cells_by="cluster",
  ordering_type = 'cluster_row_col',
  flip_percentage_mean = FALSE,
  axis_order = 'group_marker'
){
cds = cluster_cells(cds,  k = opt_k)

word_ids = as.data.frame(fData(cds)) %>%
  tibble::rownames_to_column() %>% 
  filter(rowname %in% markers) %>%
  pull(rowname)

exprs_mat = reshape2::melt(t(as.matrix(exprs(cds)[word_ids, ])))
colnames(exprs_mat) = c('Cell', 'Gene', 'Expression')
exprs_mat$Gene = as.character(exprs_mat$Gene)

if(group_cells_by == "cluster"){
  exprs_mat$Group = cds@clusters@listData$UMAP$cluster_result$optim_res$membership
} else {
  exprs_mat$Group =  colData(cds)$category
}

exprs_mat = exprs_mat %>% dplyr::filter(is.na(Group) == FALSE)

ExpVal = exprs_mat %>% dplyr::group_by(Group, Gene) %>%
  dplyr::summarize(mean = mean(Expression),
                   percentage = sum(Expression > 0) /
                     length(Expression))
ExpVal = ExpVal %>% group_by(Group) %>% mutate(avg_mean = mean(mean), avg_perc = mean(percentage))

if(flip_percentage_mean == FALSE){
  major_axis = 1
  minor_axis = 2
} else if (flip_percentage_mean == TRUE){
  major_axis = 2
  minor_axis = 1
}

res = reshape2::dcast(ExpVal[, 1:4], Group ~ Gene,
                       value.var = colnames(ExpVal)[2 + major_axis])
group_id = res[, 1]
res = res[, -1]
row.names(res) = group_id

if(ordering_type == 'cluster_row_col') {
  row_dist = stats::as.dist((1 - stats::cor(t(res)))/2)
  row_dist[is.na(row_dist)] = 1
  
  col_dist = stats::as.dist((1 - stats::cor(res))/2)
  col_dist[is.na(col_dist)] = 1
  
  ph = pheatmap::pheatmap(res,
                           useRaster = T,
                           cluster_cols=TRUE,
                           cluster_rows=TRUE,
                           show_rownames=F,
                           show_colnames=F,
                           clustering_distance_cols=col_dist,
                           clustering_distance_rows=row_dist,
                           clustering_method = 'ward.D2',
                           silent=TRUE,
                           filename=NA)
  
  ExpVal$Gene = factor(ExpVal$Gene,
                        levels = colnames(res)[ph$tree_col$order])
  ExpVal$Group = factor(ExpVal$Group,
                         levels = row.names(res)[ph$tree_row$order])
  
} else if(ordering_type == 'maximal_on_diag'){
  
  order_mat = t(apply(res, major_axis, order))
  max_ind_vec = c()
  for(i in 1:nrow(order_mat)) {
    tmp = max(which(!(order_mat[i, ] %in% max_ind_vec)))
    max_ind_vec = c(max_ind_vec, order_mat[i, tmp])
  }
  max_ind_vec = max_ind_vec[!is.na(max_ind_vec)]
  
  if(major_axis == 1){
    max_ind_vec = c(max_ind_vec, setdiff(1:length(markers), max_ind_vec))
    ExpVal$Gene = factor(ExpVal$Gene ,
                          levels = dimnames(res)[[2]][max_ind_vec])
  }
  else{
    max_ind_vec = c(max_ind_vec, setdiff(1:length(unique(exprs_mat$Group)),
                                          max_ind_vec))
    ExpVal$Group = factor(ExpVal$Group,
                           levels = dimnames(res)[[1]][max_ind_vec])
  }
} else if(ordering_type == 'none'){
  ExpVal$Gene = factor(ExpVal$Gene, levels = markers)
}

if(flip_percentage_mean){
  g = ggplot(ExpVal, aes(y = Gene,  x = Group)) +
    geom_point(aes(colour = percentage,  size = mean)) +
    viridis::scale_color_viridis(name = 'percentage') +
    scale_size(name = 'mean', range = c(0, 10))
} else {
  g = ggplot(ExpVal, aes(y = Gene,  x = Group)) +
    geom_point(aes(colour = mean,  size = percentage)) +
    viridis::scale_color_viridis(name = 'mean') +
    scale_size(name = 'percentage', range = c(1, 8))
}

if (group_cells_by == "cluster"){
  g = g + xlab("Cluster")
} else if (group_cells_by == "partition") {
  g = g + xlab("Partition")
} else{
  g = g + xlab(group_cells_by)
}

monocle_theme_opts = function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}

g = g + ylab("Word") + monocle_theme_opts() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
if(axis_order == 'marker_group') {
 g = g + coord_flip()
}

g
}

plot_words_by_group(main, markers, opt_k = 50)


