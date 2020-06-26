rm(list=ls())
library(Matrix)
library(foreach)
library(doParallel)
library(igraph)
library(dplyr)
library(rlang) # to unquote by using "!!" in "calculate.p_values" function

#baseDir <- "/users/cdursun/rice/"
#baseDir <- "C:/Users/bigjp/OneDrive/Desktop/TEST_PGR_FILES/"
baseDir <- "C:/Users/bigjp/OneDrive/Documents/PhenoGeneRankerPackage/PhenoGeneRanker/CagatayMethods/"
#testDir <- "github_test/"
testDir <- ""

inputDir <- paste0(baseDir, "datasets/networks/", testDir)
outputDir <- paste0(baseDir, "datasets/output/", testDir)

setwd(baseDir)
#source("R/github/PhenoGeneRankerFunctions.R")
#source("C:/Users/bigjp/OneDrive/Desktop/TEST_PGR_FILES/PhenoGeneRankerFunctions.R")
source("C:/Users/bigjp/OneDrive/Documents/PhenoGeneRankerPackage/PhenoGeneRanker/R/PhenoGeneRankerFunctions.R")
#source("R/github/Utils.R")
#source("C:/Users/bigjp/OneDrive/Desktop/TEST_PGR_FILES/Utils.R")

seedFile <- "tSeedGenes.txt"

# Walk matrix ID to be used
WM_ID <- 2

# Number of Random Seed Vector that will be use for p-value calculation
S <- 10
num.cores <- detectCores()
t <- Sys.time()
funcs <- as.vector(lsf.str())

# Load WM and Connectivity
loaded_vars <- load(file = paste0(inputDir,"WM_", WM_ID, ".rda"))

cat("Running for Walk Matrix: ", WM_ID, " - ", WM_Layers, "\n")

#dont need
#CandidateGenesExist <- CandidateGenes[CandidateGenes %in% gene_pool_nodes_sorted]
#cat("Number of Candidate Genes Exist : ", length(CandidateGenesExist), "\n\n")

t1 <- Sys.time()
# Read Seed Genes Used for RWR
SeedList <- read.seeds(paste0(inputDir, seedFile), gene_pool_nodes_sorted, cult_pool_nodes_sorted)
GeneSeeds <- SeedList[[1]]
CultSeeds <- SeedList[[2]]

Seeds_Score <- get.seed.scores(GeneSeeds,CultSeeds, Parameters$eta, LG, LC, Parameters$tau/LG, Parameters$phi/LC)

if (sum(Seeds_Score$Score)!=1) stop("ERROR: Seeds Scores don't add up to 1!")

# Run RWR
Random_Walk_Results <- Random_Walk_Restarts(WM, Parameters$r, Seeds_Score)
RWGeneRankDF <- rank_proteins(N, LG, Random_Walk_Results, GeneSeeds)
RWCultRankDF <- rank_cultivars(N, LG, M, LC, Random_Walk_Results, CultSeeds)


# Generate random seeds
RandomSeeds <- generate.random.seed.vector(outputDir, WM_ID, GeneSeeds, GeneConnectivity,
                                           CultSeeds, CultConnectivity, S, no.groups.gene = 2, no.groups.cult = 2)
#RandomSeeds
# calculate random ranks
Rand_Seed_Gene_Rank <- Random_Walk_Restarts_Batch(WM, RandomSeeds[["gene"]], RandomSeeds[["cult"]],
                                                  N, LG, LC,
                                                  Parameters$eta,Parameters$tau/LG, Parameters$phi/LC, Parameters$r,
                                                  funcs, 12)

# calculate p-values using random ranks
dfRank <- calculate.p_values(RWGeneRankDF, Rand_Seed_Gene_Rank, outputDir, num.cores)

# categoraiz ethe results for  candidate genes
dfRank <- categorize.ranked.genes(CandidateGenes, dfRank)


# cropped version will have gene, rank, type, p100 columns
dfRankCropped <- dfRank[, c(1:4, (ncol(dfRank)-2) )]

write.table(dfRankCropped, paste0(outputDir,"Result.txt"),sep="\t",col.names = T,
            row.names = FALSE, dec=".",quote=FALSE)




