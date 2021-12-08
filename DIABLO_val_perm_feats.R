# In this R script we create 2 data sets for validation of the method Diablo

setwd("~/Documents/Int2/Project")

library(reshape2)
library(tidyverse)

###  ### functions ### ###

# this function creates vector w (loadings) 
# it takes 3 numbers. X is a number for the seed to generate the vector, K is the number is variables in the vector, s is the persent of zeros in the vector 
# i.e. (33, 100, 0.2), 33 is the seed, 100 observations , 0.2 zeros 
# and returns this vector
weight_creation <- function(x, k, s) {
  set.seed(x)
  z <- 1-s
  w<- sample (c(0,1), size=k, replace=T, prob=c(s,z)) ### ### HERE TO MALIPULATE SPARSITY
  
  # then we change the 1 values to random from -1 to 1
  rows_with_1 = (1:length(w))[w==1] # find the values with 1
  set.seed(x)
  options(digits=2)
  w[rows_with_1] = runif(n = length(rows_with_1), min = -0.3:-0.2, max = 0.2:0.3) # replace the 1 values with random between -1 and 1
  return(w)
}

#################################################################
### ### ### ### select general variables

n=24   # samples
k=100  # variables
e=0.28
s=0
c=0

###### select from these values
# noise
#e = 0 # for No noise
#e = 0.17 # for 20% noise
#e = 0.212 # for 50% noise
#e = 0.255 # for 100% noise

# sparsity 
#s = 0 
#s = 0.2  
#s = 0.4 
#s = 0.8 

# Correlation of U vectors
#c= 0   #for cor = 1
#c= 5  #for cor = 0.8
#c= 10 #for cor = 0.5
#c= 20   #for cor = 0.2

######### select seeds
if (s==0) {wseeds <- c(2,21,5,46)
}else if(s==0.2) {wseeds <- c(2,23,5,34)
}else if(s==0.4) {wseeds <-c(2,4,5,7)
}else {wseeds <-c(88,4,82,7)}

if (c==0) {sseeds <-c(1,14,2,90,85,99) 
} else if(c==5) {sseeds <-c(1,14,2,281,485,9) 
}else if(c==10) {sseeds <-  c(1,14,2,140,384,5)
}else {sseeds <- c(1,14,2,170,485,13) 
}


###################################################################
############################# Create the data set


####### u11
set.seed(sseeds[1]) 
u11a <- rnorm(n/2, mean=-2, sd=2) # 24 numbers between 0 and 5
set.seed(sseeds[2]) 
u11b <- rnorm(n/2, mean=2, sd=2) # 24 numbers between 0 and 5
u11 <- c(u11a,u11b)
u11 <- u11 - mean(u11) # center vector

### ###  u21
set.seed(sseeds[3]) # s0.4 = 2
u21error <- runif(n = n, min = 0, max = c) ### ### ### HERE TO MANIPULATE U21 ERROR, MANIP THE COR BETWEEN U11 AND U21
u21<- u11 + u21error
u21 <- u21 - mean(u21) # center vector
#cor(u11,u21) ## Check u1 and u2 vector correlation

############ w11
w11 <- weight_creation(wseeds[1], k, s) # s0.4 =5 # s0.8 =13
w11 <- w11/sqrt(sum((w11)^2)) # scale it to length 1

######## w21
w21 <- weight_creation(wseeds[3], k, s) # s0.4 =2 # s0.8 =16
w21 <- w21/sqrt(sum((w21)^2)) # scale it to length 1

###################
X1 <- outer(u11, w11)
X2 <- outer(u21, w21) 
#hist(X1, breaks = 50,main = "Count Distribution") 
#hist(X2, breaks = 50,main = "Count Distribution") 

##########
set.seed(2)
E1 <- matrix(rnorm(n*k, mean=0, sd=e),nrow=24)### ### HERE TO MALIPULATE X1 ERROR
set.seed(3)
E2 <- matrix(rnorm(n*k, mean=0, sd=e),nrow=24) ### ### HERE TO MALIPULATE X2 ERROR

# Matrix error added checks
### Sum of Squared Variance for X1
#X1 <- as.data.frame(X1)
#X2 <- as.data.frame(X2)
#E1 <- as.data.frame(E1)
#E2 <- as.data.frame(E2)
# percentage of X1 added error
#sum(as.numeric(lapply(E1,var))^2) /sum(as.numeric(lapply(X1,var))^2)
# X2 error
#sum(as.numeric(lapply(E2,var))^2) / sum(as.numeric(lapply(X2,var))^2)

# add the residual metrices
X1 <- X1 + E1
X2 <- X2 + E2

colnames(X1) <- paste("trans", seq(1, 100) , sep="")
colnames(X2) <-  paste( "metab", seq(1, 100), sep="")


#DIABLO
library(mixOmics)

X<- list(transcriptomics = X1, metabolomics = X2)
#Y <- c("wt","pd","wt","pd","wt","pd","wt","pd","wt","pd","wt","pd" ,"wt","pd","wt","pd","wt","pd","wt","pd","wt","pd","wt","pd")
#Y <- c("wt","wt","wt","wt","wt","wt","wt","wt","wt","wt","wt","wt" ,"pd","pd","pd","pd","pd","pd","pd","pd","pd","pd","pd","pd")
Y <- c(-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
Y <- as.factor(Y)

list.keepX <- list(transcriptomics = sum(w11>0 | w11<0), metabolomics = sum(w21>0 | w21<0))

############# change design matrix DIABLO NULL
#a<- c(0,1)
#design <- as.data.frame(a)
#design$b <-  c(1,0)

########### run DIABLO

MyResult.diablo <- block.splsda(X, Y, scale = FALSE,  ncomp = 1, keepX=list.keepX)#, design = design)
options(digits=20)

permute_col <- function(col) {
  new_col <- sample(col,replace=FALSE )
  return(new_col) }

#######################################################################################################################################################
#################################### permutations X1
permlist <- list()
permlist[[1]] <- X1
for (i in 1:100){ 
  set.seed(i)
  X1perm <- apply(X1, 2, permute_col)
  permlist[[1+i]] <- X1perm
}

permres <- as.data.frame(1)

for (i in 1:100){ 
  X1 <-  permlist[[i]]
  X<- list(transcriptomics = X1, metabolomics = X2)
  
  list.keepX <- list(transcriptomics = c(sum(w11>0 | w11<0)), metabolomics = c(sum(w21>0 | w21<0)))
  
  MyResult.diablo <- block.splsda(X, Y, scale = FALSE,  ncomp = 1, keepX=list.keepX)#, design = design)
  
  corMat <- circosPlot(MyResult.diablo, cutoff=0.9999, size.labels = 1,  size.variables = 0.8)
  temp<- melt(abs(corMat))
  temp2 <- temp[order(temp$value, decreasing = TRUE),]
  temp3 <- temp2[substr(temp2$Var1, start = 1, stop = 5) != substr(temp2$Var2, start = 1, stop = 5),]
  temp3 <- head(temp3, 20)
  top20cor <- temp3$value
  permres[1:20,i] <- top20cor
  print(i)
}
colnames(permres) <- c("Top20cor")


#######################################################################################################################################################
#################################### permutataions X2
permlist <- list()
permlist[[1]] <- X2
for (i in 1:100){ 
  set.seed(i)
  X2perm <- apply(X2, 2, permute_col)
  permlist[[1+i]] <- X2perm
}

permres <- as.data.frame(1)

for (i in 1:100){ 
  X2 <-  permlist[[i]]
  X<- list(transcriptomics = X1, metabolomics = X2)
  
  list.keepX <- list(transcriptomics = c(sum(w11>0 | w11<0)), metabolomics = c(sum(w21>0 | w21<0)))
  
  MyResult.diablo <- block.splsda(X, Y, scale = FALSE,  ncomp = 1, keepX=list.keepX)#, design = design)
  
  corMat <- circosPlot(MyResult.diablo, cutoff=0.9999, size.labels = 1,  size.variables = 0.8)
  temp<- melt(abs(corMat))
  temp2 <- temp[order(temp$value, decreasing = TRUE),]
  temp3 <- temp2[substr(temp2$Var1, start = 1, stop = 5) != substr(temp2$Var2, start = 1, stop = 5),]
  temp3 <- head(temp3, 20)
  top20cor <- temp3$value
  permres[1:20,i] <- top20cor
  print(i)
}
#colnames(permres) <- c("Top20cor")

######################################################################################################################################################################
############## permutations Y

permdf <- as.data.frame(Y)
for (i in 1:100){ 
  set.seed(1+i)
  Yperm <-sample(Y, replace=FALSE, size=24)
  permdf[,1+i] <- Yperm
}
#sum(colSums(permdf) ==12)
permres <- as.data.frame(1)

for (i in 1:100){ 
  Y <- permdf[,i]
  Y <- as.factor(Y)
  
  list.keepX <- list(transcriptomics = c(sum(w11>0 | w11<0)), metabolomics = c(sum(w21>0 | w21<0)))
  
  MyResult.diablo <- block.splsda(X, Y, scale = FALSE,  ncomp = 1, keepX=list.keepX)#, design = design)
  
  corMat <- circosPlot(MyResult.diablo, cutoff=0.9999, size.labels = 1,  size.variables = 0.8)
  temp<- melt(abs(corMat))
  temp2 <- temp[order(temp$value, decreasing = TRUE),]
  temp3 <- temp2[substr(temp2$Var1, start = 1, stop = 5) != substr(temp2$Var2, start = 1, stop = 5),]
  temp3 <- head(temp3, 20)
  top20cor <- temp3$value
  permres[1:20,i] <- top20cor
  print(i)
}
#colnames(permres) <- c("Top20cor")


########################################################################



############### results
sum((abs(u11) -abs(MyResult.diablo[["variates"]][["transcriptomics"]][,1]))^2) +  
  sum((abs(u21) - abs(MyResult.diablo[["variates"]][["metabolomics"]][,1]))^2) 

sum((abs(w11) - abs(MyResult.diablo[["loadings"]][["transcriptomics"]][,1]))^2) + 
  sum((abs(w21) - abs(MyResult.diablo[["loadings"]][["metabolomics"]][,1]))^2) 

############## permutations Y
#Y <- c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1)
# create Y
#uy <- c(-0.707,-0.707,-0.707,-0.707,-0.707,-0.707,-0.707,-0.707,-0.707,-0.707,-0.707,-0.707,0.707,0.707,0.707,0.707,0.707,0.707,0.707,0.707,0.707,0.707,0.707,0.707)
#wy <- c(-0.707, 0.707)
#Y<- outer(uy, wy) 

permdf <- as.data.frame(seq(1:24), nrow=1)
for (i in 1:1000){ 
  set.seed(1+i)
  Yperm <-sample(Y, replace=FALSE, size=24)
  permdf[,i] <- Yperm
}
sum(colSums(permdf) ==12)
permres <- as.data.frame(1)
options(digits=20)

for (i in 1:1000){ 
  Y <- permdf[,i]
  Y <- as.factor(Y)
  
  list.keepX <- list(transcriptomics = c(sum(w11>0 | w11<0)), metabolomics = c(sum(w21>0 | w21<0)))
  
  MyResult.diablo <- block.splsda(X, Y, scale = TRUE,  ncomp = 1, keepX=list.keepX)
  
  serr <- sum((abs(u11) -abs(MyResult.diablo[["variates"]][["transcriptomics"]][,1]))^2) + 
    sum((abs(u21) - abs(MyResult.diablo[["variates"]][["metabolomics"]][,1]))^2) 
  
  lerr <- sum((abs(w11) - abs(MyResult.diablo[["loadings"]][["transcriptomics"]][,1]))^2) + 
    sum((abs(w21) - abs(MyResult.diablo[["loadings"]][["metabolomics"]][,1]))^2) 
  
  cerr <- cor(MyResult.diablo[["variates"]][["transcriptomics"]], MyResult.diablo[["variates"]][["metabolomics"]]) + cor(MyResult.diablo[["variates"]][["transcriptomics"]], MyResult.diablo[["variates"]][["Y"]]) + cor(MyResult.diablo[["variates"]][["metabolomics"]], MyResult.diablo[["variates"]][["Y"]])
  
  permres[i,1] <- serr
  permres[i,2] <- lerr
  permres[i,3] <- cerr
}
colnames(permres) <- c("Escore", "Eload", "Ecor")
plot(permres$Escore)

max(permres$Ecor)

# for score
hist(permres$Escore, breaks = 20, main = "Score error Distribution of permutated Y for 50% noise", xlim = c(2,8), xlab = "Score error") + 
  abline(v=3.03, col="red", lwd=3, lty=2)

# for laoad
hist(permres$Eload, breaks = 20, main = "Loading error Distribution of permutated Y for 50% noise", xlim = c(0.6,1.5), xlab = "Score error") + 
  abline(v=0.69, col="red", lwd=3, lty=2)

##################
# those with score less than 2.4 correct Y
# have loaginds
hist(permres$Eload[permres$Escore < 3], breaks = 10, main = "Loading error of min scores", xlim = c(0.5,1.5), xlab = "Loading error") + 
  abline(v=0.69, col="red", lwd=3, lty=2)

hist(permres$Escore[permres$Eload < 0.69], breaks = 10, main = "Score error of min loadings", xlim = c(2.2,5), xlab = "Score error") + 
  abline(v=3, col="red", lwd=3, lty=2)

which(permres$Escore < 3)

a<- permdf[,339]
a<- replace(a, a==0, 2)
plot(u11, col = a) +  abline(v=12.5, col="grey", lwd=1) +title("Correct vs Low error Y vector")


############### resutls ans plots
Noise <- c(0, 0.2, 0.5, 1)
Ser <- c(0, 1.5, 2.4, 3.4)
Ler <- c( 0, 0.55, 0.83, 1.2)
plot(Noise,Ser, col="red", type="l", lwd=3, )
lines(Noise,Ler,col="green", lwd=3)

a<- as.data.frame(cbind(Noise, Ser, Ler))
v<- c("scores","scores","scores","scores","loadings","loadings","loadings","loadings")
a<- as.data.frame(cbind(rbind(cbind(Noise, Ser), cbind(Noise, Ler)),v)  )
colnames(a) <- c("Noise", "er", "Error")

ggplot(a, aes(x=Noise)) + 
  geom_line(aes(y = Ser), color = "darkred", lwd=1.5) + 
  geom_line(aes(y = Ler), color="steelblue", lwd=1.5) + ggtitle("Error per noise") +
  xlab("Noise") + ylab("Error") 


############### individual checks

### ### Check scores
scores<- as.data.frame(cbind(MyResult.diablo[["variates"]][["transcriptomics"]], MyResult.diablo[["variates"]][["metabolomics"]] ))
colnames(scores) <- c("tpc1","mpc1")
scoredf <- cbind(scores, u11, u21)
head(scoredf)

par(mfrow=c(1,2))
plot(scoredf$tpc1, scoredf$u11, main ="X1 dataset", xlab = "pc1 scores (output)", ylab = "u11 (input)")
plot(scoredf$mpc1, scoredf$u21, main ="X2 dataset", xlab = "pc1 scores (output)", ylab = "u21 (input)")


scoredf$tpc1/ scoredf$u11
scoredf$mpc1/ scoredf$u21

### ### ### ### check Scores and Loadings
loadings<- as.data.frame(cbind(MyResult.diablo[["loadings"]][["transcriptomics"]], MyResult.diablo[["loadings"]][["metabolomics"]] ))
colnames(loadings) <- c("tw1","mw1")

loaddf <-as.data.frame(cbind(loadings, w11,  w21))
head(loaddf)

par(mfrow=c(1,2))
plot(loadings$tw1, w11, main = "X1 dataset", xlab = "pc1 loadings (output)", ylab = "w11 (input)")
plot(loadings$mw1, w21, main = "X2 dataset", xlab = "pc1 loadings (output)", ylab = "w21 (input)")

sum(loadings$tw1*loadings$tw1)
sum(w11*w11)
sum(loadings$tw1*w11)



###############
a <- as.data.frame(t(cbind(loaddf[,1], loaddf[,3])))
a<-as.data.frame(t(a[, colSums(a != 0) > 0]))

b <- as.data.frame(a[,1])
b$group <- "input"
colnames(b)<- c("lala", "group")
c <- as.data.frame(a[,2])
c$group <- "output"
colnames(c)<- c("lala", "group")
d <- rbind(b,c)

ggplot(d, aes(x=lala)) +  geom_density(aes(y=(..count..),group=group,color=group), adjust = 0.4)+
  geom_point(stat = "count", colour = "black") + 
  labs(title ="Loading vectors dataset 1 comp 1", x = "values of loading vectors", y = "density of obsernations") #+  ylim(50, 170)

########################################
# calculate error
#The score error for dataset 1 component 1
sum((scoredf$u11 - scoredf$tpc1)^2) # we add them because they have opposite directions
#The score error for dataset 1 component 2
sum((scoredf$u12 - scoredf$tpc2)^2) 
#The score error for dataset 2 component 1
sum((scoredf$u21 - scoredf$mpc1)^2) 
#The score error for dataset 2 component 2
sum((scoredf$u22 - scoredf$mpc2)^2) 

#total score error
sum((scoredf$u11 - scoredf$tpc1)^2) + sum((scoredf$u12 - scoredf$tpc2)^2) + sum((scoredf$u21 - scoredf$mpc1)^2) + sum((scoredf$u22 - scoredf$mpc2)^2) 

#The loading error for dataset 1 component 1
sum((loaddf$w11 - loaddf$tw1)^2)# we add them because they have opposite directions
#The score error for dataset 1 component 2
sum((loaddf$w12 - loaddf$tw2)^2)
#The score error for dataset 2 component 1
sum((loaddf$w21 - loaddf$mw1)^2)
#The score error for dataset 2 component 2
sum((loaddf$w22 - loaddf$mw2)^2) 

#total loading error
sum((loaddf$w11 - loaddf$tw1)^2) + sum((loaddf$w12 - loaddf$tw2)^2) + sum((loaddf$w21 - loaddf$mw1)^2) + sum((loaddf$w22 - loaddf$mw2)^2 )


######################################
# keep variables =/0 and rerun DIABLO

X1 <- X1[,rowSums(loadings[,1:2]) !=0]
X2 <- X2[,rowSums(loadings[,3:4]) !=0]

################################################3 Check loading order change
head(loaddf)
library(tidyverse)

w11or <- loaddf[order(loaddf$w11),] # order for in
w11or$order <- seq(1:100)
w11or <- w11or[order(w11or$tw1),] # order for out
top<- slice(w11or, 1:10)
sum(top$order <11)
bot<- slice_tail(w11or, n=10)
sum(bot$order >90)

w21or <- loaddf[order(loaddf$w21),] # order for in
w21or$order <- seq(1:100)
w21or <- w21or[order(w21or$mw1),] # order for out
top<- slice(w21or, 1:10)
sum(top$order <11)
bot<- slice_tail(w21or, n=10)
sum(bot$order >90)

############## Check cor of vectors
#################### get scores and Y
scores<- as.data.frame(cbind(MyResult.diablo[["variates"]][["transcriptomics"]], MyResult.diablo[["variates"]][["metabolomics"]], MyResult.diablo[["variates"]][["Y"]] ))
colnames(scores) <- c("t1","m1" ,"y1")

cor(scores$t1, scores$m1) + cor(scores$t1, scores$y1) + cor(scores$m1, scores$y1)

cor(MyResult.diablo[["variates"]][["transcriptomics"]], MyResult.diablo[["variates"]][["metabolomics"]]) + cor(MyResult.diablo[["variates"]][["transcriptomics"]], MyResult.diablo[["variates"]][["Y"]]) + cor(MyResult.diablo[["variates"]][["metabolomics"]], MyResult.diablo[["variates"]][["Y"]])

## plot
x1perm <- permres
x2perm <- permres
yperm <- permres

x1permnu <- permres
x2permnu <- permres
ypermnu <- permres

a<- as.data.frame(permres)
a$pairs <- seq(1:20)

df <- melt(a ,  id.vars = 'pairs', variable.name = 'series')

df2<- df %>% group_by(series) %>% summarise(mean(value))
colnames(df2) <- c("series", "mea")

plot(df2  , ylab = "correlation mean", xlab = " X2 permutations",  main = "Top 20 X1X2 feature correlations")

df3<- df2 %>% filter(mea > 0.89)
df2<- df[df$series %in% df3$series,]




