# Import packages that'll be used
library(igraph)
library(ggplot2)
library(network)
library(GGally)
library(stringr)
library(dplyr)
library(stringr)

gene_disease <- read.csv("Disease-Gene_v2015.csv", header = T, as.is = T)
miRNA_gene <- read.csv("net_miRNA_mRNA_edgelist.csv", header = F, as.is = T)
colnames(miRNA_gene) <- c("miRNA", "Gene")
colnames(gene_disease) <- c("ID", "Disease", "Type", "Gene")
gene_disease=gene_disease[gene_disease$Type!="DEG",]


#list of cellular components
genes <- unique(gene_disease$Gene)
diseases <- unique(gene_disease$Disease)
mgenes <- unique(miRNA_gene$Gene)
miRNA <- unique(miRNA_gene$miRNA)
overlap = intersect(genes, mgenes) #select overlapping genes

## Extract the rows that overlap in each dataset ##
disease_overlap= unique(gene_disease[gene_disease$Gene %in% overlap,c(2,4)]) #selected rows
ndisease = unique(gene_disease[gene_disease$Gene %in% overlap,2]) #diseases in tripartite network
miRNA_overlap = miRNA_gene[miRNA_gene$Gene %in% overlap,] #selected rows
nmiRNA=unique(miRNA_gene[miRNA_gene$Gene %in% overlap,1]) # vector of miRNAs in tripartite

## Join two data sets based on the common genes in each data set ##
tridf = inner_join(miRNA_overlap, disease_overlap) # data frame with all three parts and edges
colnames(tridf) = c("miRNA", "Gene", "Disease")

# Make two separate networks, assign all the edges a value, then union the networks together
net1=simplify(graph_from_data_frame(tridf[,1:2], directed=TRUE), remove.loops = TRUE, remove.multiple = TRUE)
net2 = simplify(graph_from_data_frame(tridf[,2:3], directed=TRUE), remove.loops = TRUE, remove.multiple = TRUE) 
V(net2)$type <- V(net2)$name %in% tridf[,2] #genes are TRUE type
V(net2)[V(net2)$type==F]

#Group synonymous diseases
for(k in 1:length(V(net2)[V(net2)$type==F]))
{
  if(k<length(V(net2)[V(net2)$type==F]))
  {
    d=V(net2)[V(net2)$type==F][k]
    n=neighbors(net2, d, mode='all')$name
    for(l in length(V(net2)[V(net2)$type==F]):1)
    {
      a=neighbors(net2, V(net2)[V(net2)$type==F][l], mode='all')$name
      if(k!=l && length(a)==length(n))
      {
        c=(a==n)
        if(all(c)==T)
        {
          newname = paste(V(net2)[V(net2)$type==F][k]$name,
                          V(net2)[V(net2)$type==F][l]$name, sep="/")
          V(net2)[V(net2)$name==d$name]$name=newname
          net2=delete_vertices(net2, V(net2)[V(net2)$type==F][l])
        }
      }
      
    }
  }
}

write.csv(V(net2)[V(net2)$type==F]$name, 'nodes.csv')
df=data.frame(head_of(net2, E(net2))$name, tail_of(net2, E(net2))$name)
write.csv(df, "cleaned_genedisease_edgelist_noDEG.csv")

