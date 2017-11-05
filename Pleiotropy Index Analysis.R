library(igraph)
library(magrittr)
library(dplyr)
library(GGally)
library(ggplot2)


###### Make two networks, experimental and predicted  #####
hmdd = readr::read_csv("HMDD raw data.csv")
hmdd = dplyr::select(hmdd, miRNA, Disease)

df = readr::read_csv("MDN_synonymsgrouped_ICDclassified3.csv")
df = dplyr::select(df, miRNA, Disease, ZScore, ICD.9.Code, ClassA)

df$miRNA %>% unique() %>% length()
hmdd$miRNA %>% unique() %>% length()

overlap = intersect(as.character(hmdd2$miRNA), as.character(df$miRNA))
overlap %>% length()

miRNA_overlap = unique(df[df$miRNA %in% overlap,]) #selected rows
hmdd_overlap = unique(hmdd[hmdd$miRNA %in% overlap,]) #selected rows

hnet = graph_from_data_frame(miRNA_overlap, directed=F)
mnet = graph_from_data_frame(hmdd_overlap, directed=F)

V(mnet)
V(hnet)

V(mnet)$type <- V(mnet)$name %in% miRNA_overlap$miRNA #miRNA are TRUE type
V(hnet)$type <- V(hnet)$name %in% hmdd_overlap$miRNA #miRNA are TRUE type


###### Plot degree distributions to compare two networks #####

mverts = V(mnet)[V(mnet)$type==T]
mdd = degree_distribution(mnet, v=mverts)
mdd_degree_dist <- data.frame("K"=c(1:length(mdd)), "Pk"=mdd)

ggplot(mdd_degree_dist, aes(x=K,y=Pk))+
  geom_point(color="black")+#geom_line()+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  xlab("k") + 
  ylab("P(k)")+
  ggtitle("Degree Distribution of miRNAs \n in Predicted Network") + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 11))


mverts = V(hnet)[V(hnet)$type==T]
mdd = degree_distribution(hnet, v=mverts)
mdd_degree_dist <- data.frame("K"=c(1:length(mdd)), "Pk"=mdd)

ggplot(mdd_degree_dist, aes(x=K,y=Pk))+
  geom_point(color="black")+#geom_line()+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  xlab("k") + 
  ylab("P(k)")+
  ggtitle("Degree Distribution of miRNAs \n in Experimental Network") + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 11))



###### Make data frames to see if the same miRNAs are pleiotropic #####

m=V(mnet)[V(mnet)$type==T]$name
d=c()
for(i in 1:length(V(mnet)[V(mnet)$type==T]))
{
  d=c(d, degree(mnet, V(mnet)[i]))
}
df=data.frame(m, as.numeric(d))
df = arrange(df, m)

m=V(hnet)[V(hnet)$type==T]$name
d=c()
for(i in 1:length(V(hnet)[V(hnet)$type==T]))
{
  d=c(d, degree(hnet, V(hnet)[i]))
}
hdf=data.frame(m, as.numeric(d))
hdf = arrange(hdf, m)

### Vectors are organized by miRNA name, ie degree in index 1 corresponds to the same miRNA in both experimental & predicted



hvec=hdf$as.numeric.d.
mvec=df$as.numeric.d.

Pxi_pred=mvec/sum(mvec)
Hx=-sum(Pxi_pred*log2(Pxi_pred))

Pxi_exp=hvec/sum(hvec)
Hx_h=-sum(Pxi_exp*log2(Pxi_exp))

plot((Pxi_pred*log2(Pxi_pred)) - (Pxi_exp*log2(Pxi_exp)))

diff = Pxi_pred - Pxi_exp
logresid = -log10(abs(resid))
test = data.frame(rbind(diff, m))
plot(diff)
ggplot(diff)+
  geom_point(color="black")+#geom_line()+
  #scale_y_continuous(trans='log10')+
  #scale_x_continuous(trans='log10')+
  xlab("miRNA index") + 
  ylab("Difference in Pleiotropy")+
  #ggtitle("") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 11))



test = full_join(df, hdf, by = 'm')
test$as.numeric.d..x = -log10(test$as.numeric.d..x / sum(test$as.numeric.d..x) )
test$as.numeric.d..y = -log10(test$as.numeric.d..y / sum(test$as.numeric.d..y))

