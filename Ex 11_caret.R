#Setup
library(tidyverse)
library(dslabs)
library(purrr)
library(pdftools)
library(lubridate)
library(magrittr)
library(broom)
library(caret) #Knn3
library(ggplot2)
library(dplyr)
library(stringr)
library(matrixStats)
library(Matrix)
library(e1071) #Confusion Matrix

require(class)
library(ElemStatLearn)
library(devtools)
library(genefilter) #***101

library(rpart)
library(quantreg)
library(randomForest)
library("Rcpp")
library("Rborist") #Loading required package: Rcpp

library(gam) #*** require "splines" and "foreach" installed #*****102

#***101
#install.packages("devtools")
#library(devtools)
#devtools::install_bioc("genefilter")
#BiocManager::install("genefilter",version = "3.8")
#install.packages("BiocManager")
#source("https://bioconductor.org/biocLite.R")
#biocLite("genefilter")
#library(genefilter)

#***102
#install.packages("splines")
#install.packages("foreach")
#the libraries are not necessary


#----------- OPTIONS ------------
options(max.print=99999) 
options(digits = 5)


# Setup RM

# Q6 - FROM RANDOM FOREST in Classification with More than Two Classes

library(rpart)
n <- 1000
sigma <- 0.25
x <- rnorm(n, 0, 1)
y <- 0.75 * x + rnorm(n, 0, sigma)
dat <- data.frame(x = x, y = y)

library(randomForest)
fit <- randomForest(y ~ ., data = dat, nodesize = 50, maxnodes = 25)
  dat %>% 
  mutate(y_hat = predict(fit)) %>% 
  ggplot() +
  geom_point(aes(x, y)) +
  geom_step(aes(x, y_hat), col = 2)

#Q1

#Solution

set.seed(1)
library(caret)
fit <- train(y ~ ., method = "Rborist",   
             tuneGrid = data.frame(predFixed = 1, 
                                   minNode = seq(25, 100, 25)),
             data = dat)
ggplot(fit, highlight = TRUE)


#Q2

#2
dat %>% 
  mutate(y_hat = predict(fit)) %>% 
  ggplot() +
  geom_point(aes(x, y)) +
  geom_step(aes(x, y_hat), col = 2)
  

#Q3

#Gmodels Q6 

set.seed(1996)
data("tissue_gene_expression")
y <- tissue_gene_expression$y
x <- tissue_gene_expression$x
str(tissue_gene_expression)
levels(tissue_gene_expression$y)
dat <- data.frame(x = x, y = y)
head(dat)

modelLookup("rpart")   
getModelInfo("rpart")
?rpart


set.seed(1996)
fit <- train(y ~ ., method = "rpart", tuneGrid = data.frame(cp = seq(0, 0.1, 0.01)), data = dat)
ggplot(fit, highlight = TRUE)
predict(fit)

fit$bestTune
fit$finalModel

#Solution

library(caret)
library(dslabs)
set.seed(1991)
data("tissue_gene_expression")

fit <- with(tissue_gene_expression, 
            train(x, y, method = "rpart",
                  tuneGrid = data.frame(cp = seq(0, 0.1, 0.01))))

ggplot(fit) 
ggplot(fit, highlight = TRUE)

#Q4

#Study the confusion matrix for the best fitting classification tree from the exercise in Q3.
#What do you observe happening for the placenta samples?

fit <- with(tissue_gene_expression, 
            train(x, y, method = "rpart",
                  tuneGrid = data.frame(cp = 0)))

confusionMatrix(y,predict(fit))

#Q5

# Rerun the analysis you did in the exercise in Q3, but this time, 
# allow rpart to split any node by using the argument control = rpart.control(minsplit = 0). 
# Look at the confusion matrix again to determine whether the accuracy increases

#Solution

set.seed(1991)
data("tissue_gene_expression")

fit_rpart <- with(tissue_gene_expression, 
                  train(x, y, method = "rpart",
                        tuneGrid = data.frame(cp = seq(0, 0.10, 0.01)),
                        control = rpart.control(minsplit = 0)))
ggplot(fit_rpart, highlight = TRUE)
confusionMatrix(fit_rpart)


#Q6

#Solution
plot(fit_rpart$finalModel)
text(fit_rpart$finalModel)

#Q7

set.seed(1991) 
dat <- data.frame(y=tissue_gene_expression['y'],tissue_gene_expression['x']) 
dat <- data.frame(x = x, y = y)
fit_rf <- train(y ~ ., data=dat, method="rf", tuneGrid=data.frame(mtry=seq(50,200,25)), nodesize=1) 
fit_rf$bestTune
ggplot(fit_rf, highlight = TRUE)


#Solution 

set.seed(1991)
library(randomForest)
fit <- with(tissue_gene_expression, 
            train(x, y, method = "rf", 
                  nodesize = 1,
                  tuneGrid = data.frame(mtry = seq(50, 200, 25))))

ggplot(fit, highlight = TRUE)

#Q8
varImp(fit)

#Q9
#Calculate the variable importance in the Random Forest call for these seven predictors 
#and examine where they rank.

#Extracting the predictor names
fit_rpart <- with(tissue_gene_expression, 
            train(x, y, method = "rpart",
                  tuneGrid = data.frame(cp = seq(0, 0.1, 0.01))))


tree_terms <- as.character(unique(fit_rpart$finalModel$frame$var[!(fit_rpart$finalModel$frame$var == "<leaf>")]))
tree_terms

z <- varImp(fit_rpart)


  getModelInfo("knn")
  modelLookup("knn")
  modelLookup("gamLoess")
  