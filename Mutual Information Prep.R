library(igraph)
library(dplyr)
library(readr)
library(stringr)
library(magrittr)
library(GGally)

##################################################################

df=read_csv("MDN_synonymsgrouped_ICDclassified3.csv")
df=dplyr::select(df, miRNA, Disease, ZScore, ICD.9.Code, ClassA)
net=graph_from_data_frame(df, directed=F)
net=simplify(net, edge.attr.comb=list(ZScore='first', ClassA='first', ICD.9.Code='first', 'ignore'))

V(net)[1:length(unique(df$miRNA))]$type <- T #miRNA are TRUE type
V(net)[(length(unique(df$miRNA))+1):length(V(net))]$type <- F 


gdf=read_csv("cleaned_genedisease_edgelist_noDEG.csv")
gnet=graph_from_data_frame(gdf[c(2:3)], directed=F)

V(gnet)[1:length(unique(gdf$Gene))]$type <- T #Genes are TRUE type
V(gnet)[(length(unique(gdf$Gene))+1):length(V(gnet))]$type <- F 

################### Thresholding ################################
net2=delete_edges(net, E(net)[E(net)$ZScore<7.7]) #Threshold till edge densities are same

edge_density(gnet)/edge_density(net2) #threshold till edge densities are same

#useless crap
#x=V(gnet)[V(gnet)$type==F]$name
#x2=V(net)[V(net)$type==F]$name
#x3=intersect(x,x2)
#write.graph(net2, "miRNA_Disease.graphml", format="graphml")
#write.graph(gnet, "gene_Disease.graphml", format="graphml")


################### Projections ################################
ddm=bipartite.projection(net2)$proj1
ddg=bipartite.projection(gnet)$proj1

#assign projections disease classes
t=data.frame(V(ddm)$name)
colnames(t)='Disease'
d=plyr::join(t,select(df, Disease, ClassA), match='first')
d=unique(d)
nrow(d)
V(ddm)$class=d$ClassA

t=data.frame(V(ddg)$name)
colnames(t)='Disease'
d=plyr::join(t,select(df, Disease, ClassA), match='first')
d=unique(d)
nrow(d)
V(ddg)$class=d$ClassA

ddm=delete_vertices(ddm, V(ddm)[degree(ddm, V(ddm))<1])
ddg=delete_vertices(ddg, V(ddg)[degree(ddg, V(ddg))<1])

#length(V(ddm)[is.na(V(ddm)$class)])/length(V(ddm))
ddm=delete_vertices(ddm, V(ddm)[is.na(V(ddm)$class)])
ddg=delete_vertices(ddg, V(ddg)[is.na(V(ddg)$class)])


write.graph(ddm, "miRNA_Disease.graphml", format="graphml")
write.graph(ddg, "gene_Disease.graphml", format="graphml")


################### Modularity and Disease class analysis ################################
gephi = read_csv("mirna_dproj_ModularityandClassdata_Gephi2.csv")
gephi = arrange(gephi, ModClass)
nrow = unique(gephi$ModClass) %>% length()
ncol = unique(gephi$ClassA) %>% length()
ntot = unique(gephi$Disease) %>% length()
classes=unique(gephi$ClassA)


tab=table(gephi)
totalvec=c()
  for(i in 0:22)
  {
    minidf = gephi[gephi$ModClass==i,] #pick one mod class
    colvec=c(rep(0,17)) #fill the first col with zeros
      for(j in 1:length(classes))
      {
        for(l in 1:nrow(minidf))
        {
          if(minidf$ClassA[l]==classes[j])
          {
            colvec[j]=colvec[j]+1
          }
        }
      }
     totalvec=c(totalvec, colvec)
  }
length(totalvec)
mat=matrix(totalvec, ncol=23)
sum(mat)
nrow(mat)
ncol(mat)
?rowsum
rowSums(mat)
row.names(mat) = classes
csums=colSums(mat)

probmat=mat
for(i in 1:ncol(mat))
{
 probmat[,i] = mat[,i]/csums[i]
}
#divide by col sums to make each value the proportion of a disease type in each module

heatmap(probmat)



