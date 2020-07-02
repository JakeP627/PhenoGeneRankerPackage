library(igraph)
library(Matrix)
library(foreach)
library(doParallel)
library(dplyr)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
#baseDir <- "/users/cdursun/rice/"
#baseDir <- "C:/Users/bigjp/OneDrive/Desktop/TEST_PGR_FILES/"
#baseDir <- "C:/Users/bigjp/Desktop/TEST_PGR_FILES/"
baseDir <- "C:/Users/bigjp/OneDrive/Documents/PhenoGeneRankerPackage/PhenoGeneRanker/CagatayMethods/"

#testDir <- "github_test/"
testDir <- ""

setwd(baseDir)
#source("R/github/PhenoGeneRankerFunctions.R")
#source("C:/Users/bigjp/OneDrive/Desktop/PhenoGeneRanker/PhenoGeneRankerFunctions.R")
#source("C:/Users/bigjp/OneDrive/Desktop/TEST_PGR_FILES/PhenoGeneRankerFunctions.R")
source("C:/Users/bigjp/OneDrive/Documents/PhenoGeneRankerPackage/PhenoGeneRanker/R/PhenoGeneRankerFunctions.R")
#source("R/github/Utils.R")
#source("C:/Users/bigjp/Desktop/TEST_PGR_FILES/Utils.R")

#inputDir <- paste0("datasets/networks/", testDir)
inputDir <- paste0("datasets/networks/", testDir)
outputDir <- paste0("datasets/output/", testDir)
setwd(paste0(baseDir, inputDir))
inputFile <- "input_files_het_par_1.txt"


t <- Sys.time()

numCores <- detectCores()
#numCores <- 12
numCores

tau <- c(1,1,1)
phi <- c(1,1,1)
myParams <- addParameters(.7,.5,.5,tau,phi,.5,.5,1,3,3)
myParams
#myWM <- createWalkMatrixUpdated(inputFile, numCores)
myWM <- createWalkMatrix(inputFile, myParams, numCores)
#myWM


length(myWM)
