
# In this R script we create 2 data sets for validation of the method JIVE

setwd("~/Documents/Int2/Project")

library(reshape2)
library(tidyverse)
library(r.jive)
library(mixOmics) # for pca

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

###############################################

n=24   # samples
k=100  # variables
e=0.04
s=0

# Noise
# for 0 = 0
# 20 =  0.015
# 50 = 0.03
# 70 =  0.046
# 80 = 0.06


##################################

##### Joint variation
# make common scores 
####### u11
set.seed(3) 
S01a <- rnorm(n/2, mean=-2, sd=2) # 24 numbers between 0 and 5
set.seed(4) 
S01b <- rnorm(n/2, mean=2, sd=2) # 24 numbers between 0 and 5
S01 <- c(S01a,S01b)
S01 <- S01 - mean(S01) # center vector
S01 <- S01/sqrt(sum((S01)^2)) # scale it to length 1

# make different loadings
######## u11
u11 <- weight_creation(4, k/2, s) # s0.4 =2 # s0.8 =16
u11 <- u11 - mean(u11)
u11 <- u11/sqrt(sum((u11)^2)) # scale it to length 1
u11 <- c(u11,rep(0,k/2))

######## u21
u21 <- weight_creation(3, k/2, s) # s0.4 =2 # s0.8 =16
u21 <- u21 - mean(u21)
u21 <- u21/sqrt(sum((u21)^2)) # scale it to length 1
u21 <- c(u21,rep(0,k/2))


############## Individual variation

####### S11
set.seed(2) 
S11 <- rnorm(n, mean=0, sd=2) # 24 numbers between 0 and 5
S11 <- S11 - mean(S11) # center vector
S11 <- S11[-n] ## make S11 orthogonal to S01 by changing the last value in the vector
S11[n] <- -sum(S01[-n] *S11[-n])/S01[n]
S11 <- S11/sqrt(sum((S11)^2)) # scale it to length 1

############ w11
w11 <- weight_creation(2, k, s) # s0.4 =5 # s0.8 =13
w11 <- w11 - mean(w11)
w11 <- w11/sqrt(sum((w11)^2)) # scale it to length 1

####### S21
set.seed(5) 
S21 <- rnorm(n, mean=0, sd=2) # 24 numbers between 0 and 5
S21 <- S21 - mean(S21) # center vector
S21 <- S21[-n] ## make S21 orthogonal to S01 by changing the last value in the vector
S21[n] <- -sum(S01[-n] *S21[-n])/S01[n]
S21 <- S21/sqrt(sum((S21)^2)) # scale it to length 1

######## w21
w21 <- weight_creation(1, k, s) # s0.4 =2 # s0.8 =16
w21 <- w21 - mean(w21)
w21 <- w21/sqrt(sum((w21)^2)) # scale it to length 1


### add everything
JointX1 <- outer(u11,S01)  * 1 # for 0.1 Joint in ttal 5% Joint var
JointX2 <- outer(u21,S01)  * 1
IndiX1 <- outer(w11,S11)   * 3
IndiX2 <-  outer(w21,S21)  * 1

X1 <- JointX1 +IndiX1 
X2 <- JointX2 + IndiX2 

############### Noise
set.seed(2)
E1 <- t(matrix(rnorm(n*k, mean=0, sd=0.046),nrow=24))### ### HERE TO MALIPULATE X1 ERROR
set.seed(3)
E2 <- t(matrix(rnorm(k*n, mean=0, sd=0.06),nrow=24)) ### ### HERE TO MALIPULATE X2 ERROR

X1 <- X1 + E1
X2 <- X2 + E2

### Check added Noise = Sum of Squared Variance for X1
X1 <- as.data.frame(t(X1))
X2 <- as.data.frame(t(X2))
E1 <- as.data.frame(t(E1))
E2 <- as.data.frame(t(E2))
# percentage of X1 added error and X2
sum(as.numeric(lapply(E1,var))) /sum(as.numeric(lapply(X1,var)))
sum(as.numeric(lapply(E2,var))) / sum(as.numeric(lapply(X2,var)))


JointX1 <- as.data.frame(t(JointX1))
JointX2 <- as.data.frame(t(JointX2))
IndiX1 <- as.data.frame(t(IndiX1))
IndiX2 <- as.data.frame(t(IndiX2))
sum(as.numeric(lapply(JointX1,var))) /sum(as.numeric(lapply(X1,var)))
sum(as.numeric(lapply(IndiX1,var))) /sum(as.numeric(lapply(X1,var)))

sum(as.numeric(lapply(JointX2,var))) /sum(as.numeric(lapply(X2,var)))
sum(as.numeric(lapply(IndiX2,var))) /sum(as.numeric(lapply(X2,var)))

# add outcome Y as X3
X3 <- as.data.frame(c(rep(-1,12), rep(1,12) ) )
X3$lala <-  c(rep(1,12), rep(-1,12) ) 
colnames(X3) <- c("var1, var2")
X3 <- t(X3)
#X3 <- as.data.frame(X3)
#X3<- X3/ sum(as.numeric(lapply(X3,var)))
#X3 <- as.matrix(X3)

Data <- list(transcriptomics  = X1, metabolomics = X2)#, outcome = X3)

##################################################
################################### Run JIVE
Results <- jive(Data,center = FALSE, scale = FALSE)

summary(Results)
showVarExplained(Results)
showHeatmaps(Results)

pcajoint1 <- pca(t(Results[["individual"]][[1]]), ncomp = 2) 
pcajoint2 <- pca(t(Results[["individual"]][[2]]), ncomp = 2) 



plot(S11,  pcajoint1[["variates"]][["X"]][,1] ,  main= "X1 scores", ylab = "JIVE Individual scores pc1", xlab = "Individual variation scores")
plot(S21,  pcajoint2[["variates"]][["X"]][,1] ,  main= "X2 scores", ylab = "JIVE Individual scores pc1", xlab = "Individual variation scores")

plot(u11, pcajoint1[["loadings"]][["X"]][,2] ,  main= "X1 loadings", ylab = "JIVE Individual loadings pc2", xlab = "Joint variation loadings")
plot(u21, pcajoint2[["loadings"]][["X"]][,2] ,  main= "X2 loadings", ylab = "JIVE Individual loadings pc2", xlab = "Joint variation loadings")

cor(u11, pcajoint1[["loadings"]][["X"]][,1] )
cor(u21, pcajoint2[["loadings"]][["X"]][,1] )

### calculate error
sum((abs(S01) -abs(as.numeric( pcajoint1[["variates"]][["X"]] )))^2) +
  sum((abs(S01) -abs(as.numeric( pcajoint2[["variates"]][["X"]] )))^2) 

sum((abs(u11) -abs(as.numeric( pcajoint1[["loadings"]][["X"]] )))^2) +
  sum((abs(u11) -abs(as.numeric( pcajoint2[["loadings"]][["X"]] )))^2) 

a<-  as.data.frame(S01)
a$or <-   pcajoint1[["variates"]][["X"]]
a$perm <- pcajoint2[["variates"]][["X"]] 
plot(a$S01, a$or, main="In vs out score",
     xlab="Input S01", ylab="Output S01 of X1") 
############## # results and plots
Noise <- c(0,20,50,70, 80, 100)
Ser <- c(0,0.0045,0.024,0.076, 1000, 1000)*10*1.7
Ler <- c(0,0.02,0.071, 0.15, 1000, 1000)*1.7
#plot(Cor,Ser, col="red", type="l", lwd=3, ylim = c(0,10))#, ylab("Error")) #+ xlab("Correlation of U1i and U2i")
#lines(Cor,Ler,col="green", lwd=3) #+ title("Erro for Changes in Cor(u1i,u2i)")# + ylab("Error") #+ xlab("Correlation of U1i and U2i")


a<- as.data.frame(cbind(Noise, Ser, Ler))
#v<- c("scores","scores","scores","scores","loadings","loadings","loadings","loadings")
#a<- as.data.frame(cbind(rbind(cbind(Cor, Ser), cbind(Cor, Ler)),v)  )
#colnames(a) <- c("Noise", "er", "Error")

ggplot(a, aes(x=Noise)) + 
  geom_line(aes(y = Ser), color = "darkred", lwd=1.5) + 
  geom_line(aes(y = Ler), color="steelblue", lwd=1.5) + ggtitle("Error per added Noise") +
  xlab("Noise percentage") + ylab("Error")  + ylim(0, 2) + xlim(0, 70)
