library(stringr)
library(magrittr)
library(GGally)
library(dplyr)


df=readr::read_csv("MDN_nodeg.csv")


DCtoICD = readr::read_csv("DiseaseConnect_to_ICD.csv") #Diseaseconnect->ICD file
DCtoICD=DCtoICD[,-1]
colnames(DCtoICD) = c("map_source", "UMLS.Code", "Disease", "ICD.9.Code", "ICD 10 Code", 
                      "Mesh Code", "OMIM Code")
DCtoICD$ICD.9.Code
ICDtoclass$ICD.9.Code

ICDtoclass = readr::read_csv(file = "ICD9_broadestlevelofhierarchy.csv") #ICD->hierarchy file
ICDtoclass=ICDtoclass[,-1]
colnames(ICDtoclass) = c("ICD.9.Code", "Name", "ClassA")
ICDtoclass$ICD.9.Code %<>% as.character()

DiseaseConnectToICDClass=dplyr::left_join(DCtoICD[3:4], ICDtoclass, copy=T)
DiseaseConnectToICDClass=DiseaseConnectToICDClass[,c("Disease", "ICD.9.Code", "ClassA")]

classifydf=dplyr::left_join(df, DiseaseConnectToICDClass) #pass 1 of classification
table(is.na(classifydf$ICD.9.Code)==T) #16256 diseases left unclassified

unclassifieddf=classifydf[is.na(classifydf$ICD.9.Code)==T,]

#need to match ICD9 codes for the combined diseases

res=c()
nums=c()
for(i in 1:nrow(unclassifieddf))
{
  strlist=data.frame(str_split(unclassifieddf[i,"Disease"], "/"))
  colnames(strlist)=c("Disease")
  strlist=dplyr::inner_join(strlist, DiseaseConnectToICDClass, by='Disease')
  v=strlist[,"ClassA"]
  c=strlist[,"ICD.9.Code"]
  newclass=v[1]
  newnum=c[1]
  res=c(res, newclass)    
  nums=c(nums, newnum)
}
?inner_join
unclassifieddf=cbind(unclassifieddf[c(2:3)], res, nums)
colnames(unclassifieddf)=c('miRNA', 'Disease', 'ClassA', "ICD9")
head(unclassifieddf)

#classifydf2= left_join(classifydf, unclassifieddf, by=c('miRNA', 'Disease'),copy=T)

unclassifieddf2=unclassifieddf
classifydf3=classifydf

for(j in 1:nrow(classifydf3))
{
  if(is.na(classifydf3$ClassA[j])==T)
  {
    xx=as.vector(unique(unclassifieddf2$ClassA[unclassifieddf2$Disease==classifydf3$Disease[j]]))
    classifydf3$ClassA[j]=xx[1]
  }
}
classifydf3$ClassA

write.csv(classifydf3, "MDN_synonymsgrouped_ICDclassified3.csv")

#write.csv(classifydf2, "MDN_synonymsgrouped_ICDclassified2.csv")
