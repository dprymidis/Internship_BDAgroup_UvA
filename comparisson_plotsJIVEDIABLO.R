
library(mixOmics) # for pca


### just data ranks 0,5,3
or<- Results
showVarExplained(or)

ort <- pca(t(or[["individual"]][[1]]), ncomp = 5) 
orm <- pca(t(or[["individual"]][[2]]), ncomp = 3) 
ort[["variates"]][["X"]][,1] 
orm[["variates"]][["X"]][,1]
ort[["loadings"]][["X"]][,1]
orm[["loadings"]][["X"]][,1]


###### X3 ranks 1,6,3,1
#resx3 <- Results
showVarExplained(resx3)

pcajx3t <- pca(t(resx3[["joint"]][[1]]), ncomp = 1) 
pcajx3m <- pca(t(resx3[["joint"]][[2]]), ncomp = 1) 
pcajx3t[["variates"]][["X"]][,1] 
pcajx3m[["variates"]][["X"]][,1]
pcajx3t[["loadings"]][["X"]][,1]
pcajx3m[["loadings"]][["X"]][,1]

##### X3 and scaled ranks 1,5,6,1
#resxc <- Results
showVarExplained(resxc)

pcaxct <- pca(t(resxc[["joint"]][[1]]), ncomp = 1) 
pcaxcm <- pca(t(resxc[["joint"]][[2]]), ncomp = 1) 
pcaxct[["variates"]][["X"]][,1]
pcaxcm[["variates"]][["X"]][,1]
pcaxct[["loadings"]][["X"]][,1]
pcaxcm[["loadings"]][["X"]][,1]

#### scale ranks 1,6,6
#ressc <- Results
showVarExplained(ressc)

pcajsct <- pca(t(ressc[["joint"]][[1]]), ncomp = 1) 
pcajscm <- pca(t(ressc[["joint"]][[2]]), ncomp = 1) 
pcajsct[["variates"]][["X"]][,1]
pcajscm[["variates"]][["X"]][,1]
pcajsct[["loadings"]][["X"]][,1]
pcajscm[["loadings"]][["X"]][,1]

# DIABLO full no keepX
MyResult.diablo

MyResult.diablo[["variates"]][["transcriptomics"]][,1] # pc1
MyResult.diablo[["variates"]][["transcriptomics"]][,2] 
MyResult.diablo[["variates"]][["metabolomics"]][,1] # pc1
MyResult.diablo[["variates"]][["metabolomics"]][,1] 

MyResult.diablo[["loadings"]][["transcriptomics"]][,1] # pc1
MyResult.diablo[["loadings"]][["transcriptomics"]][,2] 
MyResult.diablo[["loadings"]][["metabolomics"]][,1] # pc1
MyResult.diablo[["loadings"]][["metabolomics"]][,1] 


#### make all together data frames
scoresTranscriptomics <- as.data.frame(pcajx3t[["variates"]][["X"]][,1] )
colnames(scoresTranscriptomics) <- "X3"
scoresTranscriptomics$scaled <- pcajsct[["variates"]][["X"]][,1]
scoresTranscriptomics$X3scaled <- pcaxct[["variates"]][["X"]][,1]
scoresTranscriptomics$DIABLOpc1 <- MyResult.diablo[["variates"]][["transcriptomics"]][,1]
scoresTranscriptomics$DIABLOpc2 <- MyResult.diablo[["variates"]][["transcriptomics"]][,2]

plot(scoresTranscriptomics, main = "transcriptomics scores")

loadingsTranscriptomics <- as.data.frame(pcajx3t[["loadings"]][["X"]][,1] )
colnames(loadingsTranscriptomics) <- "X3"
loadingsTranscriptomics$scaled <- pcajsct[["loadings"]][["X"]][,1]
loadingsTranscriptomics$X3scaled <- pcaxct[["loadings"]][["X"]][,1]
loadingsTranscriptomics$DIABLOpc1 <- MyResult.diablo[["loadings"]][["transcriptomics"]][,1]
loadingsTranscriptomics$DIABLOpc2 <- MyResult.diablo[["loadings"]][["transcriptomics"]][,2]

plot(loadingsTranscriptomics, main = "transcriptomics loadings")

scoresMetabolomics <- as.data.frame(pcajx3m[["variates"]][["X"]][,1] )
colnames(scoresMetabolomics) <- "X3"
scoresMetabolomics$scaled <- pcajscm[["variates"]][["X"]][,1]
scoresMetabolomics$X3scaled <- pcaxcm[["variates"]][["X"]][,1]
scoresMetabolomics$DIABLOpc1 <- MyResult.diablo[["variates"]][["metabolomics"]][,1]
scoresMetabolomics$DIABLOpc2 <- MyResult.diablo[["variates"]][["metabolomics"]][,2]

plot(scoresMetabolomics, main = "Metabolomics scores")

loadingsMetabolomics <- as.data.frame(pcajx3m[["loadings"]][["X"]][,1] )
colnames(loadingsMetabolomics) <- "X3"
loadingsMetabolomics$scaled <- pcajscm[["loadings"]][["X"]][,1]
loadingsMetabolomics$X3scaled <- pcaxcm[["loadings"]][["X"]][,1]
loadingsMetabolomics$DIABLOpc1 <- MyResult.diablo[["loadings"]][["metabolomics"]][,1]
loadingsMetabolomics$DIABLOpc2 <- MyResult.diablo[["loadings"]][["metabolomics"]][,2]

plot(loadingsMetabolomics, main = "Metabolomics loadings")





## DIABLo vs X3
plot(MyResult.diablo[["variates"]][["transcriptomics"]][,1], pcajx3t[["variates"]][["X"]][,1]  ,main= "transcriptomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc1")
plot(MyResult.diablo[["variates"]][["transcriptomics"]][,2], pcajx3t[["variates"]][["X"]][,1]  ,main= "transcriptomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc2")
plot(MyResult.diablo[["variates"]][["metabolomics"]][,1], pcajx3m[["variates"]][["X"]][,1]  ,main= "metabolomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc1")
plot(MyResult.diablo[["variates"]][["metabolomics"]][,2], pcajx3m[["variates"]][["X"]][,1]  ,main= "metabolomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc2")

plot(MyResult.diablo[["loadings"]][["transcriptomics"]][,1], pcajx3t[["loadings"]][["X"]][,1]  ,main= "transcriptomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc1")
plot(MyResult.diablo[["loadings"]][["transcriptomics"]][,2], pcajx3t[["loadings"]][["X"]][,1]  ,main= "transcriptomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc2")
plot(MyResult.diablo[["loadings"]][["metabolomics"]][,1], pcajx3m[["loadings"]][["X"]][,1]  ,main= "metabolomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc1")
plot(MyResult.diablo[["loadings"]][["metabolomics"]][,2], pcajx3m[["loadings"]][["X"]][,1]  ,main= "metabolomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc2")


cor(MyResult.diablo[["variates"]][["transcriptomics"]][,1], pcajx3t[["variates"]][["X"]][,1])#  ,main= "transcriptomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc1")
cor(MyResult.diablo[["variates"]][["transcriptomics"]][,2], pcajx3t[["variates"]][["X"]][,1] )# ,main= "transcriptomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc2")
cor(MyResult.diablo[["variates"]][["metabolomics"]][,1], pcajx3m[["variates"]][["X"]][,1])#  ,main= "metabolomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc1")
cor(MyResult.diablo[["variates"]][["metabolomics"]][,2], pcajx3m[["variates"]][["X"]][,1])#  ,main= "metabolomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc2")
cor(MyResult.diablo[["loadings"]][["transcriptomics"]][,1], pcajx3t[["loadings"]][["X"]][,1])#  ,main= "transcriptomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc1")
cor(MyResult.diablo[["loadings"]][["transcriptomics"]][,2], pcajx3t[["loadings"]][["X"]][,1] )# ,main= "transcriptomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc2")
cor(MyResult.diablo[["loadings"]][["metabolomics"]][,1], pcajx3m[["loadings"]][["X"]][,1]  )#,main= "metabolomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc1")
cor(MyResult.diablo[["loadings"]][["metabolomics"]][,2], pcajx3m[["loadings"]][["X"]][,1]  )#,main= "metabolomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc2")

## DIABLo vs filtered
plot(MyResult.diablo[["variates"]][["transcriptomics"]][,1], pcajft[["variates"]][["X"]][,1]  ,main= "transcriptomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc1")
plot(MyResult.diablo[["variates"]][["transcriptomics"]][,2], pcajft[["variates"]][["X"]][,1]  ,main= "transcriptomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc2")
plot(MyResult.diablo[["variates"]][["metabolomics"]][,1], pcajfm[["variates"]][["X"]][,1]  ,main= "metabolomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc1")
plot(MyResult.diablo[["variates"]][["metabolomics"]][,2], pcajfm[["variates"]][["X"]][,1]  ,main= "metabolomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc2")


plot(MyResult.diablo[["loadings"]][["transcriptomics"]][,1], pcajft[["loadings"]][["X"]][,1]  ,main= "transcriptomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc1")
plot(MyResult.diablo[["loadings"]][["transcriptomics"]][,2], pcajft[["loadings"]][["X"]][,1]  ,main= "transcriptomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc2")
plot(MyResult.diablo[["loadings"]][["metabolomics"]][,1], pcajfm[["loadings"]][["X"]][,1]  ,main= "metabolomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc1")
plot(MyResult.diablo[["loadings"]][["metabolomics"]][,2], pcajfm[["loadings"]][["X"]][,1]  ,main= "metabolomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc2")


cor(MyResult.diablo[["variates"]][["transcriptomics"]][,1], pcajft[["variates"]][["X"]][,1])#  ,main= "transcriptomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc1")
cor(MyResult.diablo[["variates"]][["transcriptomics"]][,2], pcajft[["variates"]][["X"]][,1] )# ,main= "transcriptomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc2")
cor(MyResult.diablo[["variates"]][["metabolomics"]][,1], pcajfm[["variates"]][["X"]][,1])#  ,main= "metabolomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc1")
cor(MyResult.diablo[["variates"]][["metabolomics"]][,2], pcajfm[["variates"]][["X"]][,1])#  ,main= "metabolomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc2")
cor(MyResult.diablo[["loadings"]][["transcriptomics"]][,1], pcajft[["loadings"]][["X"]][,1])#  ,main= "transcriptomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc1")
cor(MyResult.diablo[["loadings"]][["transcriptomics"]][,2], pcajft[["loadings"]][["X"]][,1] )# ,main= "transcriptomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc2")
cor(MyResult.diablo[["loadings"]][["metabolomics"]][,1], pcajfm[["loadings"]][["X"]][,1]  )#,main= "metabolomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc1")
cor(MyResult.diablo[["loadings"]][["metabolomics"]][,2], pcajfm[["loadings"]][["X"]][,1]  )#,main= "metabolomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc2")

## DIABLo vs scale
plot(MyResult.diablo[["variates"]][["transcriptomics"]][,1], pcajsct[["variates"]][["X"]][,1]  ,main= "transcriptomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc1")
plot(MyResult.diablo[["variates"]][["transcriptomics"]][,2], pcajsct[["variates"]][["X"]][,1]  ,main= "transcriptomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc2")
plot(MyResult.diablo[["variates"]][["metabolomics"]][,1], pcajscm[["variates"]][["X"]][,1]  ,main= "metabolomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc1")
plot(MyResult.diablo[["variates"]][["metabolomics"]][,2], pcajscm[["variates"]][["X"]][,1]  ,main= "metabolomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc2")

plot(MyResult.diablo[["loadings"]][["transcriptomics"]][,1], pcajsct[["loadings"]][["X"]][,1]  ,main= "transcriptomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc1")
plot(MyResult.diablo[["loadings"]][["transcriptomics"]][,2], pcajsct[["loadings"]][["X"]][,1]  ,main= "transcriptomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc2")
plot(MyResult.diablo[["loadings"]][["metabolomics"]][,1], pcajscm[["loadings"]][["X"]][,1]  ,main= "metabolomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc1")
plot(MyResult.diablo[["loadings"]][["metabolomics"]][,2], pcajscm[["loadings"]][["X"]][,1]  ,main= "metabolomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc2")


cor(MyResult.diablo[["variates"]][["transcriptomics"]][,1], pcajsct[["variates"]][["X"]][,1])#  ,main= "transcriptomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc1")
cor(MyResult.diablo[["variates"]][["transcriptomics"]][,2], pcajsct[["variates"]][["X"]][,1] )# ,main= "transcriptomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc2")
cor(MyResult.diablo[["variates"]][["metabolomics"]][,1], pcajscm[["variates"]][["X"]][,1])#  ,main= "metabolomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc1")
cor(MyResult.diablo[["variates"]][["metabolomics"]][,2], pcajscm[["variates"]][["X"]][,1])#  ,main= "metabolomics", ylab = "JIVE Joint scores ", xlab = "DIABLO scores pc2")
cor(MyResult.diablo[["loadings"]][["transcriptomics"]][,1], pcajsct[["loadings"]][["X"]][,1])#  ,main= "transcriptomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc1")
cor(MyResult.diablo[["loadings"]][["transcriptomics"]][,2], pcajsct[["loadings"]][["X"]][,1] )# ,main= "transcriptomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc2")
cor(MyResult.diablo[["loadings"]][["metabolomics"]][,1], pcajscm[["loadings"]][["X"]][,1]  )#,main= "metabolomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc1")
cor(MyResult.diablo[["loadings"]][["metabolomics"]][,2], pcajscm[["loadings"]][["X"]][,1]  )#,main= "metabolomics", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings pc2")


###### compare 0,5,3 vs DIABLO
wherescores <- as.data.frame(ort[["variates"]][["X"]][,1]  )
colnames(wherescores) <- "pc1"
wherescores$pc2 <- ort[["variates"]][["X"]][,2] 
wherescores$pc3 <- ort[["variates"]][["X"]][,3] 
wherescores$pc4 <- ort[["variates"]][["X"]][,4] 
wherescores$pc5 <- ort[["variates"]][["X"]][,5] 
wherescores$DIABLOpc1 <- MyResult.diablo[["variates"]][["transcriptomics"]][,1]

plot(wherescores, main = "transcriptomics scores")

wherescores <- as.data.frame(ort[["loadings"]][["X"]][,1]  )
colnames(wherescores) <- "pc1"
wherescores$pc2 <- ort[["loadings"]][["X"]][,2] 
wherescores$pc3 <- ort[["loadings"]][["X"]][,3] 
wherescores$pc4 <- ort[["loadings"]][["X"]][,4] 
wherescores$pc5 <- ort[["loadings"]][["X"]][,5] 
wherescores$DIABLOpc1 <- MyResult.diablo[["loadings"]][["transcriptomics"]][,1]

plot(wherescores, main = "transcriptomics loadings")

####
wherescores <- as.data.frame(orm[["variates"]][["X"]][,1]  )
colnames(wherescores) <- "pc1"
wherescores$pc2 <- orm[["variates"]][["X"]][,2] 
wherescores$pc3 <- orm[["variates"]][["X"]][,3] 
wherescores$DIABLOpc1 <- MyResult.diablo[["variates"]][["metabolomics"]][,1]

plot(wherescores, main = "Metabolomics scores")

wherescores <- as.data.frame(orm[["loadings"]][["X"]][,1]  )
colnames(wherescores) <- "pc1"
wherescores$pc2 <- orm[["loadings"]][["X"]][,2] 
wherescores$pc3 <- orm[["loadings"]][["X"]][,3] 
wherescores$DIABLOpc1 <- MyResult.diablo[["loadings"]][["metabolomics"]][,1]

plot(wherescores, main = "Metabolomics loadings")

