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
e=0.212
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

cor(u11,u21)

#DIABLO
library(mixOmics)

X<- list(transcriptomics = X1, metabolomics = X2)
#Y <- c("wt","pd","wt","pd","wt","pd","wt","pd","wt","pd","wt","pd" ,"wt","pd","wt","pd","wt","pd","wt","pd","wt","pd","wt","pd")
#Y <- c("wt","wt","wt","wt","wt","wt","wt","wt","wt","wt","wt","wt" ,"pd","pd","pd","pd","pd","pd","pd","pd","pd","pd","pd","pd")
Y <- c(-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
cor(u21,Y)
cor(u11,Y)
Y <- as.factor(Y)

list.keepX <- list(transcriptomics = sum(w11>0 | w11<0), metabolomics = sum(w21>0 | w21<0))

############# change design matrix DIABLO NULL
a<- c(0,1)
design <- as.data.frame(a)
design$b <-  c(1,0)


MyResult.diablo <- block.splsda(X, Y, scale = FALSE,  ncomp = 1, keepX=list.keepX)#, design = design)

options(digits=20)
####################################################################################################################################################
################ results

x1x2h<- as.numeric( cor(MyResult.diablo[["variates"]][["transcriptomics"]], MyResult.diablo[["variates"]][["metabolomics"]]) )[1]
x1yh<- as.numeric( cor(MyResult.diablo[["variates"]][["transcriptomics"]], MyResult.diablo[["variates"]][["Y"]])  )[1]
x2yh<- as.numeric( cor(MyResult.diablo[["variates"]][["metabolomics"]], MyResult.diablo[["variates"]][["Y"]]) )[1]



#######################################################################################################################################################
#################################### permutations X1
permlist <- list()
permlist[[1]] <- X1
for (i in 1:100){ 
  set.seed(i)
  X1perm <- X1[sample(nrow(X1),replace=FALSE ),]
  permlist[[1+i]] <- X1perm
}

permres <- as.data.frame(1)

for (i in 1:100){ 
  X1 <-  permlist[[i]]
  X<- list(transcriptomics = X1, metabolomics = X2)
  
  list.keepX <- list(transcriptomics = c(sum(w11>0 | w11<0)), metabolomics = c(sum(w21>0 | w21<0)))
  
  MyResult.diablo <- block.splsda(X, Y, scale = FALSE,  ncomp = 2, keepX=list.keepX)#, design = design)
  
  x1x2<- as.numeric(cor(MyResult.diablo[["variates"]][["transcriptomics"]], MyResult.diablo[["variates"]][["metabolomics"]]) )[1]
  x1y <- as.numeric(cor(MyResult.diablo[["variates"]][["transcriptomics"]], MyResult.diablo[["variates"]][["Y"]])  )[1]
  x2y <- as.numeric(cor(MyResult.diablo[["variates"]][["metabolomics"]], MyResult.diablo[["variates"]][["Y"]]) )[1]
  
  permres[i,1] <- x1x2
  permres[i,2] <- x1y
  permres[i,3] <- x2y
  permres[i,4] <- x1x2 + x1y + x2y
}
colnames(permres) <- c("x1x2", "x1y", "x2y", "sum")


ggplot(permres, aes(x=x1x2)) +  geom_density(aes(x=x1x2))+
  geom_density(aes(x=x1y), colour = "red")+
  # geom_point(stat = "count", colour = "black") + 
  labs(title ="Correlation of data sets X1 perm", x = "correlation", y = "density")

#######################################################################################################################################################
#################################### permutataions X2
permlist <- list()
permlist[[1]] <- X2
for (i in 1:100){ 
  set.seed(i)
  X2perm <- X2[sample(nrow(X2),replace=FALSE ),]
  permlist[[1+i]] <- X2perm
}

permres <- as.data.frame(1)

for (i in 1:100){ 
  X2 <-  permlist[[i]]
  X<- list(transcriptomics = X1, metabolomics = X2)
  
  list.keepX <- list(transcriptomics = c(sum(w11>0 | w11<0)), metabolomics = c(sum(w21>0 | w21<0)))
  
  MyResult.diablo <- block.splsda(X, Y, scale = FALSE,  ncomp = 2, keepX=list.keepX)# ,design = design)
  
  x1x2<- as.numeric(cor(MyResult.diablo[["variates"]][["transcriptomics"]], MyResult.diablo[["variates"]][["metabolomics"]]) )[1]
  x1y <- as.numeric(cor(MyResult.diablo[["variates"]][["transcriptomics"]], MyResult.diablo[["variates"]][["Y"]])  )[1]
  x2y <- as.numeric(cor(MyResult.diablo[["variates"]][["metabolomics"]], MyResult.diablo[["variates"]][["Y"]]) )[1]
  
  permres[i,1] <- x1x2
  permres[i,2] <- x1y
  permres[i,3] <- x2y
  permres[i,4] <- x1x2 + x1y + x2y
}
colnames(permres) <- c("x1x2", "x1y", "x2y", "sum")

ggplot(permres, aes(x=x1x2)) +  geom_density(aes(x=x1x2))+
  geom_density(aes(x=x2y), colour = "red")+
  # geom_point(stat = "count", colour = "black") + 
  labs(title ="Correlation of data sets X2 perm", x = "correlation", y = "density")

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
  
  MyResult.diablo <- block.splsda(X, Y, scale = FALSE,  ncomp = 2, keepX=list.keepX)#, design = design)
  
  x1x2<- as.numeric(cor(MyResult.diablo[["variates"]][["transcriptomics"]], MyResult.diablo[["variates"]][["metabolomics"]]) )[1]
  x1y <- as.numeric(cor(MyResult.diablo[["variates"]][["transcriptomics"]], MyResult.diablo[["variates"]][["Y"]])  )[1]
  x2y <- as.numeric(cor(MyResult.diablo[["variates"]][["metabolomics"]], MyResult.diablo[["variates"]][["Y"]]) )[1]
  
  permres[i,1] <- x1x2
  permres[i,2] <- x1y
  permres[i,3] <- x2y
  permres[i,4] <- x1x2 + x1y + x2y
}
colnames(permres) <- c("x1x2", "x1y", "x2y", "sum")

ggplot(permres, aes(x=x1y)) +  geom_density(aes(x=x1y))+
  geom_density(aes(x=x2y), colour = "red")+
  # geom_point(stat = "count", colour = "black") + 
  labs(title ="Correlation of data sets Y perm", x = "correlation", y = "density")


########################################################################
############### results
sum((abs(u11) -abs(MyResult.diablo[["variates"]][["transcriptomics"]][,1]))^2) +  
  sum((abs(u21) - abs(MyResult.diablo[["variates"]][["metabolomics"]][,1]))^2) 

sum((abs(w11) - abs(MyResult.diablo[["loadings"]][["transcriptomics"]][,1]))^2) + 
  sum((abs(w21) - abs(MyResult.diablo[["loadings"]][["metabolomics"]][,1]))^2) 


## plots

library(cowplot)

########### save the values from the results 
e5x1pf <- permres # i.e. e5x1pf , error 50%, X1 data set permutation, full design
e5x2pf <- permres
e5ypf  <- permres

e5x1pn <- permres
e5x2pn <- permres
e5ypn  <- permres


a<- e5x1pf
b<- e5x1pn
c<- e5x1ph
a$group <- "full"
b$group <- "null"
c$group <- "half"
d<- rbind(a,b)

a<- e5x2pf
b<- e5x2pn
c<- e5x2ph
a$group <- "full"
b$group <- "null"
c$group <- "half"
e<- rbind(a,b)

a<- e5ypf
b<- e5ypn
c<- e5yph
a$group <- "full"
b$group <- "null"
c$group <- "half"
f<- rbind(a,b)

ggplot(d, aes(x=x1x2)) +  geom_density(alpha=.2) 

ggplot(permres, aes(x=x1x2)) +  geom_density(aes(x=x1x2))+
  geom_density(aes(x=x1y), colour = "red")

options(digits=5)

p1 <- ggplot(d, aes(x=x1x2)) +  geom_density(alpha=.2) #xlim(0, 1)  + scale_y_continuous( name = "percent")
p2<- ggplot(d, aes(x=x1y)) +  geom_density(alpha=.2) +xlim(0, 1) + scale_y_continuous( name = "percent")
p3<- ggplot(d, aes(x=x2y)) +  geom_density(alpha=.2) +xlim(0, 1)  + scale_y_continuous( name = "percent")

p4<- ggplot(e, aes(x=x1x2)) +  geom_density(alpha=.2) +xlim(0, 1) + scale_y_continuous( name = "percent")
p5<- ggplot(e, aes(x=x1y)) +  geom_density(alpha=.2) +xlim(0, 1) + scale_y_continuous( name = "percent")
p6<- ggplot(e, aes(x=x2y)) +  geom_density(alpha=.2) +xlim(0, 1) + scale_y_continuous( name = "percent")

p7<- ggplot(f, aes(x=x1x2)) +  geom_density(alpha=.2) +xlim(0, 1) + scale_y_continuous( name = "percent")
p8 <- ggplot(f, aes(x=x1y)) +  geom_density(alpha=.2) +xlim(0, 1) + scale_y_continuous( name = "percent")
p9 <- ggplot(f, aes(x=x2y)) +  geom_density(alpha=.2) +xlim(0, 1) + scale_y_continuous( name = "percent")

library(cowplot)
plot_grid(p1, p2,p3,p4,p5,p6,p7,p8,p9, ncol = 3, nrow = 3, labels=c("A",",B" ,"C", "D", "E", "F", "G", "H", "I"))

options(digits=5)
library(scales) 
ggplot(d, aes(x=x1x2,fill = group)) + geom_histogram(position = "identity", alpha = 0.2)  +xlim(0, 1)

pn1<- ggplot(d, aes(x=x1x2,fill = group)) +  geom_histogram(position = "identity", alpha = 0.2) + theme_classic() + theme(legend.position = "none") 
pn2<- ggplot(d, aes(x=x1y, fill = group)) +  geom_histogram(position = "identity", alpha = 0.2) + theme_classic()+ theme(legend.position = "none") 
pn3<- ggplot(d, aes(x=x2y, fill = group)) +  geom_histogram(position = "identity", alpha = 0.2)  + xlim(0.978, 0.983) + theme_classic()+ theme(legend.position = "none") 

pn4<- ggplot(e, aes(x=x1x2,fill = group)) +  geom_histogram(position = "identity", alpha = 0.2) + theme_classic()+ theme(legend.position = "none") 
pn5<- ggplot(e, aes(x=x1y, fill = group)) +  geom_histogram(position = "identity", alpha = 0.2)  + xlim(0.980, 0.985) + theme_classic()+ theme(legend.position = "none") # +ylim(0, 500) + theme_classic()
pn6<- ggplot(e, aes(x=x2y, fill = group)) +  geom_histogram(position = "identity", alpha = 0.2) + theme_classic()+ theme(legend.position = "none") 

pn7<- ggplot(f, aes(x=x1x2,fill = group)) +  geom_histogram(position = "identity", alpha = 0.2) + theme_classic() + xlim(0.3, 1.1)+ theme(legend.position = "none") # +ylim(0, 500) + theme_classic()
pn8 <- ggplot(f, aes(x=x1y,fill = group)) +  geom_histogram(position = "identity", alpha = 0.2) + theme_classic()+ theme(legend.position = "none") 
pn9 <- ggplot(f, aes(x=x2y,fill = group)) +  geom_histogram(position = "identity", alpha = 0.2) + theme_classic()+ theme(legend.position = "none") 


plot_grid(pn1,pn2,pn3,pn4,pn5,pn6,pn7,pn8,pn9, ncol = 3, nrow = 3,labels=c("A",",B" ,"C", "D", "E", "F", "G", "H", "I"))

########################## VIOLIN PLOTS
p<- ggplot(d, aes(x=group, y=x1x2, color= group)) + 
  geom_violin()
p + stat_summary(fun=median, geom="point", size=2, color="red")

p + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red")

##########################  BOX PLOTS
p <- ggplot(g, aes(x=group, y=x1x2)) + 
  geom_boxplot()


########################
d$perm <- "x1"
e$perm <- "x2"
f$perm <- "y"

g <- rbind(d,e,f)

save(g,file="1pc_DIABLOfull.falh.null")
#load("1pc_DIABLOfull.falh.null")


pg1<- ggplot(g, aes(x=perm, y=x1x2, fill= group)) +  geom_boxplot() + theme_classic()+ theme(legend.position = "none")
pg2<- ggplot(g, aes(x=perm, y=x1y, fill= group)) +  geom_boxplot() + theme_classic()+ theme(legend.position = "none")
pg3<- ggplot(g, aes(x=perm, y=x2y, fill= group)) +  geom_boxplot() + theme_classic()
plot_grid(pg1,pg2,pg3, ncol = 3, nrow = 1)

library(ggthemes)
g <- ggplot(mpg, aes(class, cty))
g + geom_boxplot(aes(fill=factor(cyl))) + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(title="Box plot", 
       subtitle="City Mileage grouped by Class of vehicle",
       caption="Source: mpg",
       x="Class of Vehicle",
       y="City Mileage")



###############
ggplot(d, aes(x=x1x2), group = group) +  geom_density(aes(x=x1x2))+
  #geom_density(aes(x=x1y), colour = "red")+
  # geom_point(stat = "count", colour = "black") + 
  labs(title ="Correlation of data sets Y perm", x = "correlation", y = "density")


ggplot(d, aes(x=x1x2)) +  geom_density(aes(y=(..count..),group=group,color=group), adjust = 0.4)+
  #  geom_density(aes(x=x1y), colour = "red")+
  labs(title ="Loading vectors dataset 1 comp 1", x = "values of loading vectors", y = "density of obsernations") #+  ylim(50, 170)

ggplot(d, aes(x = x1x2, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.4)

hist(d$x1x2, col='red')
hist(d$x1y, col='blue', add=TRUE)
group <- d$group
ggplot(d) + geom_density(aes(x=x1x2), color = group)#+ geom_density(aes(x=x1y)) #+ geom_density(aes(x=x2y))


geom_line(aes(x=date,y=var1),color='blue') + 
  ylab('Values')+xlab('date')




