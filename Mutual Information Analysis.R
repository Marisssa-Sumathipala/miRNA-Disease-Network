
library(stringr)
library(magrittr)
library(GGally)
library(dplyr)
library(readr)

################### Modularity and Disease class analysis ################################

gephi = read_csv("mirna_dproj_ModularityandClassdata_Gephi2.csv")
gephi = arrange(gephi, ModClass)
gephi=gephi[gephi$ClassA!="A. CONGENITAL ANOMALIES (740-759)",]
nrow = unique(gephi$ModClass) %>% length()
ncol = unique(gephi$ClassA) %>% length()
ntot = unique(gephi$Disease) %>% length()
classes=unique(gephi$ClassA)


#### Make matrix ####
totalvec=c()
for(i in 0:(nrow-1))
{
  minidf = gephi[gephi$ModClass==i,] #pick one mod class
  colvec=c(rep(0,length(classes))) #fill the first col with zeros
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
#row.names(mat) = classes
mat[1:3,1:3]


#### Make uniform probability distribution matrix ####

px=rowSums(mat) %>% as.numeric()
px = px/sum(px)
px=matrix(px, nrow = 17)

py=colSums(mat) %>% as.numeric()
py=py/sum(py)
py=matrix(py, nrow = 1)

pxy = mat/sum(mat)


pxy_px_py = pxy
for(r in 1:nrow(pxy))
{
  for(c in 1:ncol(pxy))
  {
    if(pxy_px_py[r,c]!=0)
    {
      pxy_px_py[r,c] = log(pxy[r,c]/px[r,1]/py[1,c], base=2)*pxy[r,c]
    }
  }
}

I = sum(pxy_px_py)

################### Normalize Mutual Information ################################


Hx = -sum(px*log2(px))
Hy = -sum(py*log2(py))

#Redundancy
R = I/(Hx+Hy)
Rmax = (min(Hx, Hy))/(Hx+Hy)
R/Rmax*100

#Analogous to Pearson correlation coefficient
P = I/sqrt(Hx*Hy)



