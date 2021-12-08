
# In this R script we create 2 data sets for validation of the method JIVE

setwd("~/Documents/Int2/Project")

library(reshape2)
library(tidyverse)
library(r.jive)

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
e=0.01
s=0 # constant
c=0 # constant
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
u11 <- c(u11,rep(0,50))

######## u21
u21 <- weight_creation(3, k/2, s) # s0.4 =2 # s0.8 =16
u21 <- u21 - mean(u21)
u21 <- u21/sqrt(sum((u21)^2)) # scale it to length 1
u21 <- c(u21,rep(0,50))


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

############### Noise
set.seed(2)
E1 <- t(matrix(rnorm(n*k, mean=0, sd=e),nrow=24))### ### HERE TO MALIPULATE X1 ERROR
set.seed(3)
E2 <- t(matrix(rnorm(k*n, mean=0, sd=e),nrow=24)) ### ### HERE TO MALIPULATE X2 ERROR

### Check added Noise = Sum of Squared Variance for X1
#X1 <- as.data.frame(X1)
#X2 <- as.data.frame(X2)
#E1 <- as.data.frame(E1)
#E2 <- as.data.frame(E2)
# percentage of X1 added error and X2
#sum(as.numeric(lapply(E1,var))^2) /sum(as.numeric(lapply(X1,var))^2)
#sum(as.numeric(lapply(E2,var))^2) / sum(as.numeric(lapply(X2,var))^2)

### add everything
JointX1 <- outer(u11,S01)
JointX2 <- outer(u21,S01)
IndiX1 <- outer(w11,S11)
IndiX2 <-  outer(w21,S21)

X1 <- JointX1 +IndiX1 + E1
X2 <- JointX2 + IndiX2 + E2
Data <- list(transcriptomics  = X1, metabolomics = X2)

##################################################
################################### Run JIVE
Results <- jive(Data,center = FALSE, scale = FALSE)

summary(Results)
# Visualize results
showVarExplained(Results)
# showVarExplained is also called by the "jive" S3 class default plot method
#show heatmaps
showHeatmaps(Results)

#################################################
#################################### permutations Joint X1
permlist <- list()
permlist[[1]] <- X1
for (i in 1:5){ 
  set.seed(i)
  X1perm <- X1[,sample(ncol(X1),replace=FALSE )]
  permlist[[1+i]] <- X1perm
}

joint <- as.data.frame(1)
indi <- as.data.frame(1)
res<- as.data.frame(1)

for (i in 1:3){ 
  X1 <- permlist[[i]]
  Data <- list(transcriptomics  = X1, metabolomics = X2)
  
  Results <- jive(Data, center = FALSE, scale = FALSE)
  
  temp<- summary(Results)[[3]]
  joint[1:2,i] <- temp[1,]
  indi[1:2,i] <- temp[2,]
  res[1:2,i] <- temp[3,]
}
ftable <- rbind(joint, indi, res)

mjoint <- melt(joint)
mjoint$group <- "permuted"
mjoint[mjoint$variable ==1,3] <- "original"


mjoint[mjoint$value != 0,]$group <- mjoint[mjoint$value != 0,]$variable
a<- mjoint %>% filter(value != 0)

ggplot(mjoint, aes(x=value, fill = variable)) +  geom_histogram(alpha=.4)  + #geom_density(alpha=.4, adjust = 1.5) +
  labs(title ="Estimated Joint variation for permuted X1", x = "Estimated Joint variation ", y = "Counts")

