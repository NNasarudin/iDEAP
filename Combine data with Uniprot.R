##############################################
# Read data from files                       #
# Combine preprocess data with Uniprot data  #
##############################################
df <- read.csv("no_missing.csv", header = TRUE)
uniProt <- read.csv("uniprot.csv", header = TRUE)

#Assign the file as Data
sData=df
uData=uniProt

#change to matrix form
#sMatrix<-data.matrix(sData)
#uMatrix<-data.matrix(uData)

gene.amount<-nrow(uData)
turn <- c(1:gene.amount)

d=matrix(NA, gene.amount, ncol = 50) #create new matrix



for(i in turn){
  
  # --- Get the information of uniProt data ---
  tmp.uniProt<- uData$Entry[i]
  uni=(matrix(tmp.uniProt))
  tmp.gene<- uData$Gene.names[i]
  x<-unlist(strsplit(matrix(tmp.gene), "[[:blank:]]+"))
  #print(length(x))
  
  y=length(x)+1
  #print(y)
  
  d[i,1]=uni
  
  if (length(x)==1){
    d[i,2]=x
  }
  d[i,2:y]=x
  
  
}
write.csv(d,"new_uni.csv") #save into csv file

newdf <- read.csv("new_uni.csv", header = TRUE)
new=newdf

genes.amount<-nrow(sData)
turn<- c(1:genes.amount)
turns<- c(1:gene.amount)

c=matrix(NA,gene.amount, ncol=2) #create new matrix
#b=matrix(NA,genes.amount, ncol=1)
for(i in turn){
  
  tmpData<-sData$X[i]
  
  
    for(a in turns){
      if (tmpData==d[a,2:50]){
          
        dfData<-uData$Entry[a]
        pro<-(matrix(dfData))
        c[a,1]=pro
        
        gene<-new$V2[a]
        gene1<-(matrix(gene))
        c[a,2]=gene1
        
        #exp<-sData$GSM93789[a]
        #exp1<-(matrix(exp))
        #c[a,2]=exp1
      
    } 
  }
   
}

write.csv(c,"new_pro.csv") #save into csv file


nw<- read.csv("new_pro.csv", header = TRUE)
nw1<-nw


####remove missing data
r=which(nw1$V1%in%NA)
newd=(nw1[-r,])
write.csv(newd,"no_missing_pro.csv")

newd <- read.csv("no_missing_pro.csv", header = TRUE)
nd<-as.matrix(newd)
nd=nd[,3:4]
e=matrix(NA,nrow(nd), ncol=26) #create new matrix
######################################################################################
###combined uniID with gene expression
######################################################################################
turn1<-(nrow(nd))
for(i in turn){
  
    tmpData<-sData$X[i]
    #print(tmpData)
    for(a in 1:12145){
      #print("okt")
      if (tmpData==nd[a,2]){
        #print("ok")
        exp<-sData[i,]
        
        e[a,1]=nd[a,1]
        e[a,2]=nd[a,2]
        e[a,3]=sData$X[i]
        e[a,4]=sData$GSM93789[i]
        e[a,5]=sData$GSM93920[i]
        e[a,6]=sData$GSM93921[i]
        e[a,7]=sData$GSM93922[i]
        e[a,8]=sData$GSM93923[i]
        e[a,9]=sData$GSM93924[i]
        e[a,10]=sData$GSM93925[i]
        e[a,11]=sData$GSM93926[i]
        e[a,12]=sData$GSM93927[i]
        e[a,13]=sData$GSM93928[i]
        e[a,14]=sData$GSM93929[i]
        e[a,15]=sData$GSM93932[i]
        e[a,16]=sData$GSM93938[i]
        e[a,17]=sData$GSM93939[i]
        e[a,18]=sData$GSM93941[i]
        e[a,19]=sData$GSM93943[i]
        e[a,20]=sData$GSM93944[i]
        e[a,21]=sData$GSM93946[i]
        e[a,22]=sData$GSM93948[i]
        e[a,23]=sData$GSM93950[i]
        e[a,24]=sData$GSM93952[i]
        e[a,25]=sData$GSM93954[i]
      }
    }
  
}
write.csv(e,"final_data.csv")


