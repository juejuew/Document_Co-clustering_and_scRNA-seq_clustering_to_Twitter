#################################################################################
#################################################################################
#################### Seurat (Tweets from four distinct users) ###################
#################################################################################
#################################################################################

setwd("C:/Users/jueju/Desktop/LDA and Co-clustering algorithms with data")
getwd()

library(dplyr)
library(useful)
library(Seurat)
library(patchwork)
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
doc_metadata = doc_metadata %>% dplyr::select(Keywords, Doc_Text, screen_name, user_id_new)
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

######## Create Seurat Object #########
pbmc = CreateSeuratObject(
  word_expression,
  project = "CreateSeuratObject",
  assay = "Twitter",
  names.field = 1,
  meta.data = doc_metadata
)

######## Pre-processing ########

# Normalizing the data
pbmc <- NormalizeData(pbmc)

#scale the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Identification of highly variable features 
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 1321)

#Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npc = 100)

#####################################
############ scree plot #############
#####################################
Seurat_PCA = pbmc@reductions$pca@stdev
Seurat_PCA = (Seurat_PCA^2)/sum(Seurat_PCA^2)

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

############ Distance matrix ############## 

# (1) Choice 1: Ferg's distance matrix
distance_matrix <- read.csv("Robyn_tweets_distance.csv", header = FALSE)

# (2) Choice 2: Euclidean distance matrix
coord = pbmc[["pca"]]@cell.embeddings
distance_matrix  = dist(coord, diag = T, upper = T)

k_range_s = c(5:20,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,150,170,200, 300, 500,700)
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


########################################################
############# Useful Visualization Tools ###############
########################################################
k = k_range_s[which.max(avg_widths_s)]
pbmc <- FindNeighbors(pbmc, dims = 1:50, k.param = k)
pbmc <- FindClusters(pbmc, resolution = 0.5)

### UMAP ###
# Clustering results
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 5)

### Violin plot ###
VlnPlot(pbmc, features = c("storm", "health", "goblu", "vegan"))

### Feature plot ###
FeaturePlot(pbmc, features = c("storm", "health", "goblu", "vegan"))


