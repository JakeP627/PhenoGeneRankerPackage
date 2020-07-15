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
baseDir <- "C:/Users/bigjp/OneDrive/Documents/PhenoGeneRankerPackage/PhenoGeneRanker/"

#testDir <- "github_test/"
testDir <- ""

setwd(baseDir)
#source("R/github/PhenoGeneRankerFunctions.R")
#source("C:/Users/bigjp/OneDrive/Desktop/PhenoGeneRanker/PhenoGeneRankerFunctions.R")
#source("C:/Users/bigjp/OneDrive/Desktop/TEST_PGR_FILES/PhenoGeneRankerFunctions.R")
#source("C:/Users/bigjp/OneDrive/Documents/PhenoGeneRankerPackage/PhenoGeneRanker/R/PhenoGeneRankerFunctionsCorrected.R")
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

#myParams <- addParameters(.5,.5,.5,1,3,3)
#myParams

# if (r > 1 || r <= 0){ stop("Incorrect r, it must be between 0 and 1")}
# if (sum(tau)/LG != 1) {stop(sprintf("Incorrect tau, the sum of its component divided by %d must be 1", LG))}
# if (sum(phi)/LP != 1) {stop(sprintf("Incorrect phi, the sum of its component divided by %d must be 1", LP))}
# if (eta > 1 || eta < 0){ stop("Incorrect eta, it must be between 0 and 1")}

#myWM <- createWalkMatrixUpdated(inputFile, numCores)
source("C:/Users/bigjp/OneDrive/Documents/PhenoGeneRankerPackage/PhenoGeneRanker/R/PhenoGeneRankerFunctions.R")
myWM <- CreateWalkMatrix(inputFile, numCores)
RWR <- RandomWalkRestarts(myWM, c("LOC_Os06g39750","LOC_Os09g29820"),c(), TRUE, 12)
#length(RWR)
#myWM


length(myWM)
