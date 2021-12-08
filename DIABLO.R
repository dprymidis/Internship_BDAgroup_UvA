
library(readr)
library(tidyverse)
library(mixOmics)

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
mdata<-  mdata[-c(25:27),]
tdata<-  tdata[-c(25:27),]
rownames(tdata) <- rownames(mdata)


X <- list(transcriptomics = tdata, 
          metabolomics = mdata)
Y <- rownames(mdata)
Y <- c("Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl" ,"Pdef","Pdef","Pdef","Pdef","Pdef","Pdef","Pdef","Pdef","Pdef","Pdef","Pdef","Pdef")
Y <- as.factor(Y)

#Y <- sapply(Y, function(x){substring(x, 4, 4)})


#list.keepX <- list(transcriptomics = c(2000, 2000), metabolomics = c(100,100))
list.keepX <- list(metabolomics = 20, transcriptomics = 20)
#list.keepX <- list(metabolomics = c(8,827), transcriptomics = c(40,31784))

MyResult.diablo <- block.splsda(X, Y, ncomp = 1, keepX=list.keepX)#, ncomp = 1)#, scale = FALSE)
#saveRDS(MyResult.diablo, "DIABLO_result")

plotIndiv(MyResult.diablo,  legend=TRUE,legend.title = "Groups")## sample plot

#plotVar(MyResult.diablo, var.names = c(FALSE, FALSE),  legend=TRUE, pch=c(16,16))

plotDiablo(MyResult.diablo, ncomp = 1)

plotVar(MyResult.diablo, cutoff = 0.842) ## variable plot
plotLoadings(MyResult.diablo, comp = 1, contrib = "max", ndisplay=20)

X11()
circosPlot(MyResult.diablo, cutoff=0.9, size.labels = 1.2,  size.variables = 0.6)
X11()
cimDiablo(MyResult.diablo, color.blocks = c('darkorchid', 'brown1'), comp = 1, legend.position = "topright")



#################### get scores and loadings
scores<- as.data.frame(cbind(MyResult.diablo[["variates"]][["transcriptomics"]], MyResult.diablo[["variates"]][["metabolomics"]] ))
colnames(scores) <- c("tpc1", "tpc2","mpc1", "mpc2")


##### transcriptomics
trans.load<- as.data.frame(abs(MyResult.diablo[["loadings"]][["transcriptomics"]]))
trans.load$id <- rownames(trans.load)
#trans.vae.exp <- MyResult.diablo[["prop_expl_var"]][["transcriptomics"]]
#trans.load[,1] <- trans.load[,1] * trans.vae.exp[1]
#trans.load[,2] <- trans.load[,2] * trans.vae.exp[2]
#trans.load$import.loads <- trans.load[,1] + trans.load[,2]
trans.important.variables <-as.data.frame( trans.load[order(trans.load$comp1, decreasing =  TRUE),] )

ITAG4 <- read.table("Data/ITAG4.0_descriptions.txt", sep = "\t",  quote = "")
colnames(ITAG4) <- c("ID","Desc")
trans.important.variables$id <- gsub("\\..*","",rownames(trans.important.variables))
ITAG4$id <- gsub("\\..*","",ITAG4$ID)
trans.important.variables <- merge(trans.important.variables, ITAG4, by = "id")
trans.important.variables <-as.data.frame( trans.important.variables[order(trans.important.variables$comp1, decreasing =  TRUE),] )

vars<- head( trans.important.variables, 20)
vars$id<- NULL
vars$comp1 <-NULL
colnames(vars)<- c("ID", "Description")
write.table(vars, "./mydata.txt", sep="\t")

#### metabolomics
metab.load<- as.data.frame(abs(MyResult.diablo[["loadings"]][["metabolomics"]]))
metab.vae.exp <- MyResult.diablo[["prop_expl_var"]][["metabolomics"]]
metab.load[,1] <- metab.load[,1] * metab.vae.exp[1]
metab.load[,2] <- metab.load[,2] * metab.vae.exp[2]
metab.load$import.loads <- metab.load[,1] + metab.load[,2]

metab.important.variables <- metab.load[order(metab.load$import.loads, decreasing =  TRUE),]

############## Check cor of vectors
#################### get scores and Y
scores<- as.data.frame(cbind(MyResult.diablo[["variates"]][["transcriptomics"]], MyResult.diablo[["variates"]][["metabolomics"]], MyResult.diablo[["variates"]][["Y"]] ))
colnames(scores) <- c("t1", "t2","m1", "m2","y1", "y2" )

cor(scores$t1, scores$m1)
cor(scores$t1, scores$y1)
cor(scores$m1, scores$y1)

cor(scores$t2, scores$m2)
cor(scores$t2, scores$y2)
cor(scores$m2, scores$y2)


#################3 check corrleation of 

a<- trans.important.variables$comp1 !=0
a<- trans.important.variables[a,]

b<- metab.important.variables$comp1 !=0
b<- metab.important.variables[b,]

tnames <- a$id
mnames <- rownames(b)

selm <- mdata[,colnames(mdata) %in% mnames]
colnames(tdata)<- gsub("\\..*","",colnames(tdata))
selt <- tdata[,colnames(tdata) %in% tnames]

sum(colnames(tdata) %in% tnames)

selvars <- cbind(selm, selt)

plot(selvars$`195.95414 Da 73.30 s`, selvars$Solyc01g090210)
plot(selvars$`278.03968 Da 46.85 s`, selvars$Solyc11g067280)
plot(selvars$`334.06618 Da 81.02 s`, selvars$Solyc10g085850)
plot(selvars$Solyc01g090210 , selvars$Solyc06g069470)

colnames(selvars)


