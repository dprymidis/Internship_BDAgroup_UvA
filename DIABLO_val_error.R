# In this R script we create 2 data sets for validation of the method Diablo

setwd("~/Documents/Int2/Project")
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
e=0.263
s=0
c=0

###### select from these values
# noise
#e = 0 # for No noise
#e = 0.355 # for 20% noise
#e = 0.447 # for 50% noise
#e = 0.535 # for 100% noise

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

if (c==0) {sseeds <-c(1,8,2,90,85,99) 
} else if(c==5) {sseeds <-c(1,8,2,11,485,9) 
}else if(c==10) {sseeds <-  c(1,14,2,16,384,12)
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

############## u12
set.seed(sseeds[4]) # 85
u12 <- rnorm(n, mean=0, sd=2) # 24 numbers between 0 and 5
u12<- u12[-n]
#u12 <- u12 - mean(u12) # center vector
u12[n] <- -sum(u11[-n] *u12)/u11[n] # makes u12 orthogonal to u11 b adding the "correct" nth observation
u12 <- u12 - mean(u12) # center vector

############### u22 
set.seed(sseeds[6]) 
u22error <- runif(n = n, min = 0, max = c) ### ### ###HERE TO MANIPULATE U22 ERROR, MANIP THE CORR BETWEEN U12 AND U22
u22<- u12 + u22error
u22 <- u22 - mean(u22) # center it
cor(u12, u22)
u22 <- u22[-n] ## make u22 orthogonal to u21 by changing the last value in the vector
u22[n] <- -sum(u21[-n] *u22)/u21[n]
u22 <- u22 - mean(u22) # center it again

u11<- u11*2
u21<- u21*2
############ w11
w11 <- weight_creation(wseeds[1], k, s) # s0.4 =5 # s0.8 =13
w11 <- w11/sqrt(sum((w11)^2)) # scale it to length 1

############# w12
w12 <- weight_creation(wseeds[2], k-1, s) # k=100 this is 14
w12[k] <- -sum(w11[-k] *w12)/w11[k] # make it orthogonal
w12<- w12/sqrt(sum(w12^2)) # make the length 1 to be rthonormal

######## w21
w21 <- weight_creation(wseeds[3], k, s) # s0.4 =2 # s0.8 =16
w21 <- w21/sqrt(sum((w21)^2)) # scale it to length 1

############# w22
w22 <- weight_creation(wseeds[4], k-1, s) # s0.4 =14 # s0.8 =24
w22[k] <- -sum(w21[-k] *w22)/w21[k]
w22<- w22/sqrt(sum(w22^2)) # make the length 1 to finally make it orthonormal

###################

X1 <- outer(u11, w11) + outer(u12, w12) 
X2 <- outer(u21, w21) + outer(u22, w22) 
#hist(X1, breaks = 50,main = "Count Distribution") 
#hist(X2, breaks = 50,main = "Count Distribution") 

##########
set.seed(2)
E1 <- matrix(rnorm(n*k, mean=0, sd=e),nrow=24)### ### HERE TO MALIPULATE X1 ERROR
set.seed(3)
E2 <- matrix(rnorm(n*k, mean=0, sd=e),nrow=24) ### ### HERE TO MALIPULATE X2 ERROR

X1 <- X1  + E1
X2 <- X2  + E2

# Matrix error added checks
### Check added Noise = Sum of Squared Variance for X1
X1 <- as.data.frame(X1)
X2 <- as.data.frame(X2)
E1 <- as.data.frame(E1)
E2 <- as.data.frame(E2)
# percentage of X1 added error and X2
sum(as.numeric(lapply(E1,var))^2) /sum(as.numeric(lapply(X1,var))^2)
sum(as.numeric(lapply(E2,var))^2) / sum(as.numeric(lapply(X2,var))^2)



#DIABLO
library(mixOmics)

X<- list(transcriptomics = X1, metabolomics = X2)
Y <- c("wt","wt","wt","wt","wt","wt","wt","wt","wt","wt","wt","wt" ,"pd","pd","pd","pd","pd","pd","pd","pd","pd","pd","pd","pd")
# <- c(-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
Y <- as.factor(Y)

list.keepX <- list(transcriptomics = c(sum(w11>0 | w11<0),sum(w12>0 | w12<0)), metabolomics = c(sum(w21>0 | w21<0),sum(w22>0 | w22<0)))


MyResult.diablo <- block.splsda(X, Y, scale = FALSE,  ncomp = 2, keepX=list.keepX)

sum((abs(u11) -abs(MyResult.diablo[["variates"]][["transcriptomics"]][,1]))^2) + 
  sum((abs(u12) - abs(MyResult.diablo[["variates"]][["transcriptomics"]][,2]))^2) + 
  sum((abs(u21) - abs(MyResult.diablo[["variates"]][["metabolomics"]][,1]))^2) + 
  sum((abs(u22) - abs(MyResult.diablo[["variates"]][["metabolomics"]][,2]))^2) 

sum((abs(w11) - abs(MyResult.diablo[["loadings"]][["transcriptomics"]][,1]))^2) + 
  sum((abs(w12) - abs(MyResult.diablo[["loadings"]][["transcriptomics"]][,2]))^2) + 
  sum((abs(w21) - abs(MyResult.diablo[["loadings"]][["metabolomics"]][,1]))^2) + 
  sum((abs(w22) - abs(MyResult.diablo[["loadings"]][["metabolomics"]][,2]))^2)


### ### Check scores
scores<- as.data.frame(cbind(MyResult.diablo[["variates"]][["transcriptomics"]], MyResult.diablo[["variates"]][["metabolomics"]] ))
colnames(scores) <- c("tpc1", "tpc2","mpc1", "mpc2")
scoredf <- cbind(scores, u11, u12,u21, u22)
head(scoredf)

par(mfrow=c(2,2))

plot(scoredf$tpc1, scoredf$u11, main ="X1 dataset", xlab = "pc1 scores (output)", ylab = "u11 (input)")
plot(scoredf$tpc2, scoredf$u12, main ="X1 dataset", xlab = "pc2 scores (output)", ylab = "u12 (input)")
plot(scoredf$mpc1, scoredf$u21, main ="X2 dataset", xlab = "pc1 scores (output)", ylab = "u21 (input)")
plot(scoredf$mpc2, scoredf$u22, main ="X2 dataset", xlab = "pc2 scores (output)", ylab = "u22 (input)")

scoredf$tpc1/ scoredf$u11
scoredf$tpc2/ scoredf$u12
scoredf$mpc1/ scoredf$u21
scoredf$mpc2/ scoredf$u22

cor(scoredf$tpc2,u12)
### ### ### ### check Scores and Loadings
loadings<- as.data.frame(cbind(MyResult.diablo[["loadings"]][["transcriptomics"]][1:100,], MyResult.diablo[["loadings"]][["metabolomics"]][1:100,] ))
colnames(loadings) <- c("tw1", "tw2", "mw1", "mw2")

loaddf <-as.data.frame(cbind(loadings, w11, w12, w21, w22))
head(loaddf)

par(mfrow=c(2,2))
plot(loadings$tw1, w11, main = "X1 dataset", xlab = "pc1 loadings (output)", ylab = "w11 (input)")
plot(loadings$tw2, w12, main = "X1 dataset", xlab = "pc1 loadings (output)", ylab = "w12 (input)")
plot(loadings$mw1, w21, main = "X1 dataset", xlab = "pc1 loadings (output)", ylab = "w21 (input)")
plot(loadings$mw2, w22, main = "X1 dataset", xlab = "pc1 loadings (output)", ylab = "w22 (input)")

sum(loadings$tw1*loadings$tw1)
sum(w11*w11)
sum(loadings$tw1*w11)



############### compare in and out loading vectors
a <- as.data.frame(t(cbind(loaddf[,1], loaddf[,5])))
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

######################################## Calculate error individually
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
sum((loaddf$w12 + loaddf$tw2)^2)
#The score error for dataset 2 component 1
sum((loaddf$w21 - loaddf$mw1)^2)
#The score error for dataset 2 component 2
sum((loaddf$w22 + loaddf$mw2)^2) 

#total loading error
sum((loaddf$w11 - loaddf$tw1)^2) + sum((loaddf$w12 - loaddf$tw2)^2) + sum((loaddf$w21 - loaddf$mw1)^2) + sum((loaddf$w22 - loaddf$mw2)^2 )


############## # results and plots

#
Noise <- c(0,20,50,100)
Scoreer <- c(0,16,26,39)
Loadinger <- c(0,0.32,0.47, 0.63)*10
#plot(Cor,Ser, col="red", type="l", lwd=3, ylim = c(0,10))#, ylab("Error")) #+ xlab("Correlation of U1i and U2i")
#lines(Cor,Ler,col="green", lwd=3) #+ title("Erro for Changes in Cor(u1i,u2i)")# + ylab("Error") #+ xlab("Correlation of U1i and U2i")


a<- as.data.frame(cbind(Noise, Scoreer, Loadinger))
#v<- c("scores","scores","scores","scores","loadings","loadings","loadings","loadings")
#a<- as.data.frame(cbind(rbind(cbind(Cor, Ser), cbind(Cor, Ler)),v)  )
#colnames(a) <- c("Noise", "er", "Error")

ggplot(a, aes(x=Noise)) + 
  geom_line(aes(y = Ser), color = "darkred", lwd=1.5) + 
  geom_line(aes(y = Ler), color="steelblue", lwd=1.5) + ggtitle("Error per added Noise") +
  xlab("Noise percentage") + ylab("Error") # + ylim(0, 5)
  


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

w12or <- loaddf[order(loaddf$w12),] # order for in
w12or$order <- seq(1:100)
w12or <- w12or[order(w12or$tw2),] # order for out
top<- slice(w12or, 1:10)
sum(top$order <11)
bot<- slice_tail(w12or, n=10)
sum(bot$order >90)

w22or <- loaddf[order(loaddf$w22),] # order for in
w22or$order <- seq(1:100)
w22or <- w21or[order(w22or$mw2),] # order for out
top<- slice(w22or, 1:10)
sum(top$order <11)
bot<- slice_tail(w22or, n=10)
sum(bot$order >90)




###################### DIABLO add 900 random variables to each data set

X1 <- outer(u11, w11) + outer(u12, w12) 
X2 <- outer(u21, w21) + outer(u22, w22) 

# random variables
set.seed(15)
X1rv <- matrix(rnorm(n*900, mean=0, sd=1),nrow=24)
set.seed(22)
X2rv <- matrix(rnorm(n*900, mean=0, sd=1),nrow=24) 

X1 <- cbind(X1, X1rv)
X2 <- cbind(X2, X2rv)

set.seed(2)
E1rv <- matrix(rnorm(n*1000, mean=0, sd=e),nrow=24)### ### HERE TO MALIPULATE X1 ERROR
set.seed(3)
E2rv <- matrix(rnorm(n*1000, mean=0, sd=e),nrow=24) ### ### HERE TO MALIPULATE X2 ERROR
X1 <- X1 + E1rv
X2 <- X2 + E2rv

#X1 <- X1*5
#X2 <- X2 + 1

X<- list(transcriptomics = X1, metabolomics = X2)
Y <- c("wt","wt","wt","wt","wt","wt","wt","wt","wt","wt","wt","wt" ,"pd","pd","pd","pd","pd","pd","pd","pd","pd","pd","pd","pd")
Y <- as.factor(Y)


list.keepX <- list(transcriptomics = c(20,20), metabolomics = c(20,20))
library(mixOmics)


MyResult.diablo <- block.splsda(X, Y, scale = TRUE,  ncomp = 2)#, keepX=list.keepX)

plotIndiv(MyResult.diablo) ## sample plot
plotDiablo(MyResult.diablo, ncomp = 1)
plotVar(MyResult.diablo, cutoff = 0.8) ## variable plot
plotLoadings(MyResult.diablo, comp = 1, contrib = "max", size.name = 1, ndisplay = 10)
a<- as.data.frame(order(abs(w11)))
a<- as.data.frame(order(w11))
b<- as.data.frame(order(abs(w21)))

circosPlot(MyResult.diablo, cutoff=0.9)
cimDiablo(MyResult.diablo, color.blocks = c('darkorchid', 'brown1'), comp = 1)#, margin=c(8,20), legend.position = "right")


order(MyResult.diablo[["loadings"]][["transcriptomics"]][,1])
