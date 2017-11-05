# Import packages that'll be used
library(igraph)
library(ggplot2)
library(network)
library(GGally)
library(stringr)
library(dplyr)
library(stringr)

gene_disease <- readr::read_csv("cleaned_genedisease_edgelist_noDEG.csv")
miRNA_gene <- readr::read_csv("net_miRNA_mRNA_edgelist.csv")
colnames(miRNA_gene) <- c("miRNA", "Gene")





#list of cellular components
genes <- unique(gene_disease$Gene)
diseases <- unique(gene_disease$Disease)
mgenes <- unique(miRNA_gene$Gene)
miRNA <- unique(miRNA_gene$miRNA)
overlap = intersect(genes, mgenes) #select overlapping genes

## Extract the rows that overlap in each dataset ##
disease_overlap= unique(gene_disease[gene_disease$Gene %in% overlap,2:3]) #selected rows
miRNA_overlap = miRNA_gene[miRNA_gene$Gene %in% overlap,] #selected rows

## Join two data sets based on the common genes in each data set ##
tridf = inner_join(miRNA_overlap, disease_overlap) # data frame with all three parts and edges
colnames(tridf) = c("miRNA", "Gene", "Disease")

# Make two separate networks, assign all the edges a value, then union the networks together
net1=simplify(graph_from_data_frame(tridf[,1:2], directed=TRUE), remove.loops = TRUE, remove.multiple = TRUE)
net2 = simplify(graph_from_data_frame(tridf[,2:3], directed=TRUE), remove.loops = TRUE, remove.multiple = TRUE) 
V(net2)$type <- V(net2)$name %in% gene_disease[,2] #genes are TRUE type

net3=graph.union(net1, net2)
trinet = net3

# vector of all the miRNA vertices
tri_miRNAs = V(net3)[1:nrow(unique(miRNA_overlap[,1]))]

tail(tri_miRNAs)


# Find all mirNA-Disease connections
vlist_noDEG = list()
for(a in 1:length(tri_miRNAs))
{
  vert=tri_miRNAs[a]
  a_neighbors = neighbors(net3, vert, mode="out")
  cvect = c(vert)
  for(i in 1:length(a_neighbors))
  {
    c = neighbors(net3, a_neighbors[i], mode="out")
    cvect = c(cvect, c)
  }
  vlist_noDEG=append(vlist_noDEG, list(cvect))
}





# Calculate mean and standard deviation using binomial distribution

#mean = (Km*Kd)/Ngenes; Ngenes = length(overlap)
#stdev = sqrt(mean*(1-p))
#p=Kd/Ngenes
#n=Km
means_noDEG=list()
slist_noDEG=list()
numgenes = length(overlap)
for(index in 1:length(vlist_noDEG))
{
  m=vlist_noDEG[[index]][1]
  dlist = unique(vlist_noDEG[[index]][-1])
  km=degree(trinet, m)
  mus=c(km)
  stdevs=c(km)
  for(dis in 1:length(dlist))
  {
    d=dlist[dis]
    kd=degree(trinet, d)
    p = kd/numgenes
    mu=(kd*km)/numgenes
    mus = c(mus, mu)
    sd = sqrt(mu*(1-p))
    stdevs = c(stdevs,sd)
  }
  means_noDEG = append(means_noDEG, list(mus))
  slist_noDEG = append(slist_noDEG, list(stdevs))
}



# Put means and stdevs in a data frame
msdf_noDEG=data.frame()
for(e in 1:length(means_noDEG))
{
  mi = c(rep(attr(means_noDEG[[e]][1], "name"), length(means_noDEG[[e]])-1))
  dl=attr(means_noDEG[[e]][-1], "name")
  l = as.double(means_noDEG[[e]][-1])
  s=as.double(slist_noDEG[[e]][-1])
  mudf = data.frame(mi, dl, l,s)
  msdf_noDEG = rbind(msdf_noDEG, mudf)
}
names(msdf_noDEG)=c("miRNA", "Disease", "Means", "Stdev")

# Put miRNA, Disease, and Freq in a data frame
frqdf_noDEG=list()
for (val in 1:length(means_noDEG))
{
  lis=table(attr(vlist_noDEG[[val]][-1], "name"))
  mir = c(rep(attr(means_noDEG[[val]][1], "name"), length(unique(vlist_noDEG[[val]]))-1))
  minidf=data.frame(mir, lis)
  frqdf_noDEG=rbind(frqdf_noDEG, minidf)
}
colnames(frqdf_noDEG)=c("miRNA", "Disease", "Frequency")

newdf = join(frqdf_noDEG, msdf_noDEG)

zlist=c()
for(y in 1:length(newdf[,1]))
{
  zs=(newdf[y,3]-newdf[y,4])/newdf[y,5]
  zlist=c(zlist, zs)
}

newdf=cbind(newdf, zlist)
colnames(newdf)=c("miRNA", "Disease", "Freq", "Mean", "Stdev", "ZScore")
write.csv(newdf, file="MDN_nodeg.csv")

