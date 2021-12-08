

X1 <- t(X1)
X2 <- t(X2)

X<- list(transcriptomics = X1, metabolomics = X2)
Y <- c("wt","wt","wt","wt","wt","wt","wt","wt","wt","wt","wt","wt" ,"pd","pd","pd","pd","pd","pd","pd","pd","pd","pd","pd","pd")
Y <- as.factor(Y)

list.keepX <- list(transcriptomics = 50, metabolomics = 50)


MyResult.diablo <- block.splsda(X, Y, scale = FALSE,  ncomp = 1, keepX=list.keepX)

MyResult.diablo[["variates"]][["transcriptomics"]][,1]
MyResult.diablo[["variates"]][["metabolomics"]][,1]

plot(MyResult.diablo[["variates"]][["transcriptomics"]][,1], S01,main= "X1 dataset", ylab = "Real scores ", xlab = "DIABLO scores")
plot(MyResult.diablo[["variates"]][["metabolomics"]][,1], S01,main= "X2 dataset", ylab = "Real scores ", xlab = "DIABLO scores")
     
plot(MyResult.diablo[["loadings"]][["transcriptomics"]][,1], u11,main= "X1 dataset", ylab = "Real loadings ", xlab = "DIABLO loadings")
plot(MyResult.diablo[["loadings"]][["metabolomics"]][,1], u21,main= "X2 dataset", ylab = "Real loadings ", xlab = "DIABLO loadings")


plot(MyResult.diablo[["variates"]][["transcriptomics"]][,1],  pcajoint1[["variates"]][["X"]][,1] ,main= "X1 dataset", ylab = "JIVE Joint scores ", xlab = "DIABLO scores")
plot(MyResult.diablo[["variates"]][["metabolomics"]][,1],  pcajoint2[["variates"]][["X"]][,1] ,  main= "X2 dataset", ylab = "JIVE Joint scores", xlab = "DIABLO scores")

plot(MyResult.diablo[["loadings"]][["transcriptomics"]][,1],  pcajoint1[["loadings"]][["X"]][,1] ,main= "X1 dataset", ylab = "JIVE Joint loadings ", xlab = "DIABLO loadings")
plot(MyResult.diablo[["loadings"]][["metabolomics"]][,1],  pcajoint2[["loadings"]][["X"]][,1] ,  main= "X2 dataset", ylab = "JIVE Joint loadings", xlab = "DIABLO loadings")

