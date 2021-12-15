

library(readr)
library(tidyverse)
library(ggpubr)
library(reshape2)
library(matrixTests)
setwd("~/Documents/Int2/Project")


# check either metabolomics or transcriptomics data
mdata <- read_csv("Data/processed/metabolomics_polar_data_vst_new.csv")
mdata<- mdata %>% column_to_rownames(var="...1") # make first col colnames
mdata <- mdata[ , grepl( "wt" , names( mdata ) ) ] # we take only the wild type samples
mfeatures <- rownames(mdata) # save the feature names
sample.names <- colnames(mdata)
mdata <- as.data.frame(t(mdata)) # samples in rows and features in columns


################################################## check assumptions
##### Check normality
ggqqplot(mdata$`140.95096 Da 12.21 s`)  +ggtitle("Normality check qqplot metabolomics" )

shapiro.test(mdata$`123.95933 Da 12.27 s`)


#  we have to remove features with zero variance.
a<- as.vector(sapply(mdata, function(x) var(x)) == 0)
data <- mdata[,!a]

shapirotres <- as.data.frame(sapply(data, function(x) shapiro.test(x))  )
shapirotres <- as.data.frame(t(shapirotres))
sum(shapirotres$p.value > 0.05)
#  The data is normal if the p-value is above 0.05.

# check homoscedasticity
levtest<- data

levtest$treatment <- rownames(levtest)
levtest$treatment <- sapply(levtest$treatment, function(x){substring(x, 1, 6)})


melted <- melt(levtest, id.vars=daa("treatment")) # make the groups 
groups <- group_by(melted, treatment)

group<- as.vector(levtest$treatment)
levtest$treatment <- NULL


Levens.res<- col_levene(levtest, group)
# p<0.05 means unqeual variances

sum(Levens.res$pvalue > 0.05)

