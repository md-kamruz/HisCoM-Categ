rm(list=ls())
library(MASS)
library(Matrix)
library(pracma)
library(magic) 

################################################################################


source("HisCoM_Categ_source.R")

################################################################################
###############################Data#############################################
Pep <- read.table("pathway_set.txt", header = TRUE)
dat <- read.csv("Example_dat.csv", header = TRUE)

y <- dat$Phenotype
path <- Pep$group
path_var <-  Pep$var

####For cumulative logit model
res <- HisCoM_Categ_POM(y, path, path_var, data = dat, indx = NULL,
                                    maxiter = 100,  lambda1 = 8000, lambda2 = 8000,
                                    tol = 0.00001)
res$beta_coef

####For adjacent category logit model
res <- HisCoM_Categ_ACL(y, path, path_var, data = dat, indx = NULL,
                        maxiter = 100,  lambda1 = 8000, lambda2 = 8000,
                        tol = 0.00001)
res$beta_coef

####For baseline logit model
res <- HisCoM_Categ_BCL(y, path, path_var, data = dat, indx = NULL,
                        maxiter = 100,  lambda1 = 8000, lambda2 = 8000,
                        tol = 0.00001)
res$beta_coef
