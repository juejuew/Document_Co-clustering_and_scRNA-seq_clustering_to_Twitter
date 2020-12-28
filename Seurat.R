###############################################
###############################################
#################### Seurat ###################
############################################### 
############################################### 

setwd("C:/Users/jueju/Desktop/LDA and Co-clustering algorithms with data")
getwd()

library(dplyr)
library(useful)
library(dplyr)
library(Seurat)
library(patchwork)
library(cluster)

####### expression ########
sparse_matrix <- read.csv("test_new_stem.csv", header = FALSE)
sparse_matrix = as.matrix(sparse_matrix)
sparse_matrix = t(sparse_matrix)
word_expression = sparse_matrix
dim(word_expression)

####### Doc_metadata ########
doc_metadata <- read.csv("df_doc_topic_new_stemming.csv")
dim(doc_metadata)
str(doc_metadata)

doc_metadata  = doc_metadata  %>% select(Keywords, Doc_Text, screen_name, user_id_new)
dim(doc_metadata)

####### Word_metadata ########
word_metadata <- read.csv("word_metadata_new.csv")

str(word_metadata)
word_metadata = word_metadata %>% select(word) %>% mutate(gene_short_name = word)

dim(word_metadata)

######## Change row_col names#######
rownames(word_expression) <- word_metadata$word
colnames(word_expression) <- rownames(doc_metadata)

######## Create Seurat Object #########
pbmc = CreateSeuratObject(
  word_expression,
  project = "CreateSeuratObject",
  assay = "Twitter",
  names.field = 1,
  meta.data = doc_metadata
)

pbmc

# Normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)


#scale the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 1321)

#Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npc = 100)

#####################################
############ scree plot #############
#####################################
Seurat_PCA = pbmc@reductions$pca@stdev
PCA <- data.frame(
  PCA_components = 1:100,
  Seurat = Seurat_PCA
)

PCA %>% ggplot(aes(x = PCA_components, y = Seurat)) + geom_point() +
  ylab("Proportion of variance explained by each component") +
  ggtitle("Seurat")

#################################################
############# Silhouette Analysis ###############
#################################################
distance_matrix <- read.csv("Robyn_tweets_distance.csv", header = FALSE)
#NMI_Seurat3 = read.csv("NMI_Seurat3.csv")
#NMI_Seurat3 = NMI_Seurat3 %>% select(-X)

k_range_s  =c(5:20,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,150,170,200, 300, 500,700)
avg_widths_s = vector()

for (i in 1:length(k_range_s)){
  #x = as.numeric(NMI_Seurat3[,i])
  num = as.integer(k_range[i])
  pbmc <- FindNeighbors(pbmc, dims = 1:50, k.param = num)
  pbmc <- FindClusters(pbmc, resolution = 0.5)
  x = as.numeric(Idents(pbmc))
  
  si = silhouette(x, distance_matrix)
  ssi <- summary(si)
  avg_widths_s[i] = ssi$avg.width
}

max(avg_widths_s)
length(avg_widths_s)

######### Visualization #########
library(ggplot2)
summary_s = data.frame(k_range_s, as.numeric(avg_widths_s))
str(summary_s)
ggplot(data=summary_s, aes(x=k_range_s, y=avg_widths_s)) +
  geom_line()+
  geom_point() + scale_x_continuous(name="k", breaks = c(5,60,120,170,200, 300, 500,700)) +
  scale_y_continuous(name="Average Silhouette Width", limits=c(-0.3,0.05)) + 
  ggtitle("Silhouette for Seurat(num of PC = 50)") +
  theme(axis.text.x = element_text(size=8, angle=45))

