
#install.packages("r.jive")
library(r.jive)

library(readr)
library(tidyverse)

setwd("~/Documents/Int2/Project")

### ### ###  Load metabolomics data
mdata <- read_csv("Data/processed/metabolomics_polar_data_vst_new.csv")
mdata<- mdata %>% column_to_rownames(var="...1") # make first col colnames
mdata <- mdata[ , grepl( "wt" , names( mdata ) ) ] # we take only the wild type samples
mfeatures <- rownames(mdata) # save the feature names
sample.names <- colnames(mdata)
mdata <- as.data.frame(t(mdata)) # samples in rows and features in columns

### ### ### Load transcriptomic data
#the RNA seq data pre-processed (normalised ans vst transformed) -> "RNAseq_data_norm_vst.csv"
tdata <- read_csv("Data/processed/RNAseq_data_norm_vst.csv")
tdata<- tdata %>% column_to_rownames(var="...1")
tdata <- tdata[ , grepl( "WT" , names( tdata ) ) ] # we take only the wild type samples
tfeatures <- rownames(tdata)
tdata <- as.data.frame(t(tdata))

## remove rep from the data, in order to have only wild type and p deficiency 
tdata<-  tdata[-c(25:27),]
rownames(tdata) <- rownames(mdata)

# scale vars 1
mdata <- scale(mdata, center = TRUE)
tdata <- scale(tdata, center = TRUE)

mdata[is.na(mdata)] <- 0
tdata[is.na(tdata)] <- 0

metabolomics <- t(mdata)
transcriptomics <- t(tdata)

# add outcome Y as X3
outcome <- as.data.frame(c(rep(0,12), rep(1,12) ) )
outcome$lala <-  c(rep(1,12), rep(0,12) ) 
colnames(outcome) <- c("var1, var2")
outcome <- t(outcome)

Data2 <- list(transcriptomics  = transcriptomics, metabolomics = metabolomics)#, outcome = outcome)

#data(BRCA_data)

#########3 run JIVE
#Results <- jive(Data2)
Results <- jive(Data2)#,method="given",rankJ=2,rankA=c(4,2))  #, orthIndiv = FALSE)
summary(Results)

# Using BIC rank selection
#BIC_result <- jive(Data2, method="bic")
#summary(BIC_result)

# Visualize results
showVarExplained(Results)
# showVarExplained is also called by the "jive" S3 class default plot method

png("JIVE_rep.filt.png",height=665,width=805)
showHeatmaps(Results)
dev.off() 


#show PCA plots
showPCA(Results,1,c(1,1))



########################################################################
png("VarExplained.png",height=300,width=450)  
showVarExplained(Results)  
dev.off()

png("ruy.png",height=465,width=705)
showHeatmaps(Results)
dev.off() 



####################### find genes
library(mixOmics) # for pca

pcajointtrans <- pca(t(Results[["joint"]][[1]]), ncomp = 1) 
pcajointmetab <- pca(t(Results[["joint"]][[2]]), ncomp = 1) 

 pcajointtrans[["loadings"]][["X"]]
 pcajointmetab[["loadings"]][["X"]]
 
 ##### transcriptomics
 trans.load<- as.data.frame(abs(  pcajointtrans[["loadings"]][["X"]] ))
trans.load$names <- rownames(transcriptomics)
 trans.important.variables <- trans.load[order(trans.load$PC1, decreasing =  TRUE),]
 
 ITAG4 <- read.table("Data/ITAG4.0_descriptions.txt", sep = "\t",  quote = "")
 colnames(ITAG4) <- c("ID","Desc")
 trans.important.variables$id <- gsub("\\..*","",trans.important.variables$names)
 ITAG4$id <- gsub("\\..*","",ITAG4$ID)
 trans.important.variables <- merge(trans.important.variables, ITAG4, by = "id")
 
 save(trans.important.variables,file="JIVE_X3_scaled_trans.vars")
 load(JIVE_X3_trans.vars)
 #### metabolomics
 metab.load<- as.data.frame(abs(MyResult.diablo[["loadings"]][["metabolomics"]]))
 metab.vae.exp <- MyResult.diablo[["prop_expl_var"]][["metabolomics"]]
 metab.load[,1] <- metab.load[,1] * metab.vae.exp[1]
 metab.load[,2] <- metab.load[,2] * metab.vae.exp[2]
 metab.load$import.loads <- metab.load[,1] + metab.load[,2]
 
 metab.important.variables <- metab.load[order(metab.load$import.loads, decreasing =  TRUE),]
 



