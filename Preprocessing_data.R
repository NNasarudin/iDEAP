##########################################
### Change probe id to gene symbol #######
##########################################
## R 3.3.2
## Affymetrix Human Genome U133A Array

source("http://bioconductor.org/biocLite.R")
# Install packages
biocLite("affy")
biocLite("annotate")
biocLite("hgu133a.db")
# load packages
library("affy")
library("annotate")
require("org.Hs.eg.db")
library("hgu133a.db")

# hgu133a.db
#Load the file
df <- read.table("GSE4107_series_matrix.txt", header = TRUE)
#Assign the file as sDAta
sData=df
row.names(sData)=NULL
row.names(sData)=sData[,1]
sData=sData[,2:ncol(sData)]
#head(sData)

#Checking if sDAta is matric or not
is.matrix(sData)
sData=as.matrix(sData)

Entrez_IDs = unlist(mget(rownames(sData), hgu133aSYMBOL, ifnotfound=NA))
write.table(Entrez_IDs, "ID_Symbol.txt", sep = "\t")
#write.csv(Entrez_IDs,"Symbol.csv")
#Combine the symbol with matrix
data=cbind(Entrez_IDs,sData)
row.names(data) = NULL
row.names(data) = data[,1]
#Remove probe ID
data = data[,2:ncol(data)]
mode(data) <- "numeric"
write.table(data, "GDS4107_EntID.txt", sep = "\t")
write.csv(data,"GSE4107_symbol.csv",)

# to solve repeated gene IDs
name=unique(rownames(data))            # left non-repeated gene IDs
d=matrix(NA,length(name),ncol(data))   # create new matrix
colnames(d)=colnames(data)
rownames(d)=name
for(i in 1:length(rownames(d))){
  tmp=rownames(d)[i]
  id=which(rownames(data)==tmp)
  if(length(id)>1){
    tmp_matrix=matrix(NA,length(id),ncol(data))
    for(j in 1:length(id)){
      idx=id[j]
      tmp_matrix[j,]=data[idx,]
    }
    d[i,]=colMeans(tmp_matrix)   # average gene expression values
    
  }
  if(length(id)==1){
    d[i,]=data[id,]
  }
}
# print(d)  # new matrix formed
write.csv(d,"no_repeated1.csv")
# remove whole row with missing gene ID
r=which(rownames(d)%in%NA)
newd=(d[-r,])
write.csv(newd,"no_missing1.csv")








