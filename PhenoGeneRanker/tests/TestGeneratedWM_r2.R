# Modified from RWR-MH Package library(RandomWalkRestartMH)
rm(list=ls())
library(testthat)
library(Matrix)

# This is required to run using terminal @linux server before Utils.R
if (strsplit(getwd(),"/")[[1]][length(strsplit(getwd(),"/")[[1]])]=="R") setwd("..")
context("Matrices Computation")
source("R/Utils_r1.R")
source("R/RWFunctions_r2.R")



isDeleteFiles <- FALSE

if (isDeleteFiles){
  testFiles <- paste0(network_dir, c("Expr.txt", "PPI.txt", "GeneCult.txt", "EL.txt", "LTSS.txt"))
  # remove all previously created files
  do.call(file.remove, list(testFiles[do.call(file.exists, list(testFiles))]))
}

create.base.test.networks <- function(){
  
  ## Multiplex
  geneGRF1 <- igraph::graph(c("g1","g2","g1","g3","g2","g3"), directed = FALSE)
  geneDF1 <- igraph::as_data_frame(geneGRF1)
  geneDF1$weight <- 1
  saveRDS(geneDF1, file=paste0(network_dir, "Expr.rds"))
  
  geneGRF2 <- igraph::graph(c("g1","g3","g2","g3","g3","g4","g1","g4"), directed = FALSE)
  geneDF2 <- igraph::as_data_frame(geneGRF2)
  geneDF2$weight <- 1
  saveRDS(geneDF2, file=paste0(network_dir, "PPI.rds"))
  
  SupraAdjacencyMatrixExpected <- Matrix::Matrix(c(0,0.5,0.5,0,0.5,0,0,0,
                                        0.5,0,0.5,0,0,0.5,0,0,
                                        0.5,0.5,0,0,0,0,0.5,0,
                                        0,0,0,0,0,0,0,0.5,
                                        0.5,0,0,0,0,0,0.5,0.5,
                                        0,0.5,0,0,0,0,0.5,0,
                                        0,0,0.5,0,0.5,0.5,0,0.5,
                                        0,0,0,0.5,0.5,0,0.5,0),
                                      byrow = TRUE, nrow = 8, ncol = 8)
  
  colnames(SupraAdjacencyMatrixExpected) <- c("g1_1","g2_1","g3_1","g4_1","g1_2","g2_2",
                                   "g3_2","g4_2")
  rownames(SupraAdjacencyMatrixExpected) <- c("g1_1","g2_1","g3_1","g4_1","g1_2","g2_2",
                                   "g3_2","g4_2")
  
  SupraAdjacencyMatrixExpected <- as(SupraAdjacencyMatrixExpected, "dgCMatrix")
  
  
  SupraAdjacencyMatrixNormExpected <- matrix(c(0,0.3333333,0.3333333,0,0.3333333,0,0,0,
                                    0.3333333,0,0.3333333,0,0,0.5,0,0,
                                    0.3333333,0.3333333,0,0,0,0,0.25,0,
                                    0,0,0,0,0,0,0,0.3333333,
                                    0.3333333,0,0,0,0,0,0.25,0.3333333,
                                    0,0.3333333,0,0,0,0,0.25,0,
                                    0,0,0.3333333,0,0.3333333,0.5,0,0.3333333,
                                    0,0,0,1,0.3333333,0,0.25,0),
                                  byrow = TRUE, nrow = 8, ncol = 8)
  
  colnames(SupraAdjacencyMatrixNormExpected) <- c("g1_1","g2_1","g3_1","g4_1","g1_2","g2_2","g3_2","g4_2")
  rownames(SupraAdjacencyMatrixNormExpected) <- c("g1_1","g2_1","g3_1","g4_1","g1_2","g2_2","g3_2","g4_2")
  
  
  SupraAdjacencyMatrixNormExpected <- as(SupraAdjacencyMatrixNormExpected, "dgCMatrix")
  
  
  ## Multiplex-Heterogeneous
  cultGRF1 <- igraph::graph(c("c1","c3","c2","c5","c5","c4","c5","c3"), directed = FALSE)
  cultDF1 <- igraph::as_data_frame(cultGRF1)
  cultDF1$weight <- 1
  
  saveRDS(cultDF1, file=paste0(network_dir, "EL.rds"))
  

  geneCultDF <- data.frame(from=c("g1","g3"),to=c("c1","c5"))
  geneCultDF$weight <- 1
  
  saveRDS(geneCultDF, file=paste0(network_dir, "GeneStrain.rds"))
  
  
  WMExpected <- matrix(
    c(0,0.3333333,0.1666667,0,0.1666667,0,0,0,0.25,0,0,0,0,
      0.1666667,0,0.1666667,0,0,0.5,0,0,0,0,0,0,0,
      0.1666667,0.3333333,0,0,0,0,0.125,0,0,0,0,0,0.2500000,
      0,0,0,0,0,0,0,0.3333333,0,0,0,0,0,
      0.1666667,0,0,0,0,0,0.125,0.3333333,0.25,0,0,0,0,
      0,0.3333333,0,0,0,0,0.125,0,0,0,0,0,0,
      0,0,0.1666667,0,0.1666667,0.5,0,0.3333333,0,0,0,0,0.2500000,
      0,0,0,1,0.1666667,0,0.125,0,0,0,0,0,0,
      0.5000000,0,0,0,0.5000000,0,0,0,0,0,0.5,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0.1666667,
      0,0,0,0,0,0,0,0,0.50,0,0,0,0.1666667,
      0,0,0,0,0,0,0,0,0,0,0,0,0.1666667,
      0,0,0.5000000,0,0,0,0.500,0,0,1,0.5,1,0),
    byrow = TRUE, nrow = 13, ncol = 13)
  
  
  colnames(WMExpected) <- c("g1_1","g2_1","g3_1","g4_1","g1_2","g2_2","g3_2","g4_2",
                                            "c1_1","c2_1","c3_1","c4_1","c5_1")
  rownames(WMExpected) <- c("g1_1","g2_1","g3_1","g4_1","g1_2","g2_2","g3_2","g4_2",
                                            "c1_1","c2_1","c3_1","c4_1","c5_1")
  WMExpected <- as(WMExpected, "dgCMatrix")
  print("Expected matrix for BASE test is given below: ")
  printSpMatrix2(WMExpected, col.names = TRUE)  
  
  saveit( WMExpected=WMExpected, 
          SupraAdjacencyMatrixExpected=SupraAdjacencyMatrixExpected,
          SupraAdjacencyMatrixNormExpected=SupraAdjacencyMatrixNormExpected,
         file = paste0(network_dir, "WMExpected.rda"))
}
create.extended.test.networks <- function(){
  ## Multiplex
  geneGRF1 <- igraph::graph(c("g1","g2","g1","g3","g2","g3"), directed = FALSE)
  geneDF1 <- igraph::as_data_frame(geneGRF1)
  geneDF1$weight <- 1
  saveRDS(geneDF1, file=paste0(network_dir, "Expr.rds"))
  
  geneGRF2 <- igraph::graph(c("g1","g3","g2","g3","g3","g4","g1","g4"), directed = FALSE)
  geneDF2 <- igraph::as_data_frame(geneGRF2)
  geneDF2$weight <- 1
  saveRDS(geneDF2, file=paste0(network_dir, "PPI.rds"))
  
  ## Multiplex Cult
  cultGRF1 <- igraph::graph(c("c1","c2","c1","c3","c2","c3"), directed = FALSE)
  cultDF1 <- igraph::as_data_frame(cultGRF1)
  cultDF1$weight <- 1
  saveRDS(cultDF1, file=paste0(network_dir, "EL.rds"))
  
  cultGRF2 <- igraph::graph(c("c1","c3","c2","c3","c3","c4","c1","c4"), directed = FALSE)
  cultDF2 <- igraph::as_data_frame(cultGRF2)
  cultDF2$weight <- 1
  saveRDS(cultDF2, file=paste0(network_dir, "LTSS.rds"))
  
  
  geneCultDF <- data.frame(from=c("g1","g3"),to=c("c1","c5"))
  geneCultDF$weight <- 1
  

  saveRDS(geneCultDF, file=paste0(network_dir, "GeneStrain.rds"))

  WMExtExpected <- matrix(
    c(0, 0.3333333, 0.1666667, 0, 0.1666667, 0, 0, 0, 0.25, 0, 0, 0, 0, 0.2500000, 0, 0, 0, 0,
      0.1666667, 0, 0.1666667, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0.1666667, 0.3333333, 0, 0, 0, 0, 0.125, 0, 0, 0, 0, 0, 0.2500000, 0, 0, 0, 0, 0.25,
      0, 0, 0, 0, 0, 0, 0, 0.3333333, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0.1666667, 0, 0, 0, 0, 0, 0.125, 0.3333333, 0.25, 0, 0, 0, 0, 0.2500000, 0, 0, 0, 0,
      0, 0.3333333, 0, 0, 0, 0, 0.125, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0.1666667, 0, 0.1666667, 0.5, 0, 0.3333333, 0, 0, 0, 0, 0.2500000, 0, 0, 0, 0, 0.25,
      0, 0, 0, 1, 0.1666667, 0, 0.125, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0.2500000, 0, 0, 0, 0.2500000, 0, 0, 0, 0, 0, 0.3333333, 0, 0, 0.1666667, 0, 0, 0, 0, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.125, 0, 0.5, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0.125, 0, 0, 0.25, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.125, 0, 0, 0, 0.333333, 0,
      0, 0, 0.2500000, 0, 0, 0, 0.250, 0, 0, 0.5, 0.3333333, 0.5, 0, 0, 0, 0, 0, 0.50,
      0.2500000, 0, 0, 0, 0.2500000, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0.25, 0.3333333, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0.25, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3333333, 0, 0, 0.1666667, 0.5, 0, 0.3333333, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0.1666667, 0, 0.25, 0, 0,
      0, 0, 0.2500000, 0, 0, 0, 0.250, 0, 0, 0, 0, 0, 0.125, 0, 0, 0, 0, 0),
    byrow = TRUE, nrow = 18, ncol = 18)
  
  
  colnames(WMExtExpected) <- c("g1_1","g2_1","g3_1","g4_1","g1_2","g2_2","g3_2","g4_2",
                            "c1_1","c2_1","c3_1","c4_1","c5_1", "c1_2","c2_2","c3_2","c4_2","c5_2")
  rownames(WMExtExpected) <- c("g1_1","g2_1","g3_1","g4_1","g1_2","g2_2","g3_2","g4_2",
                            "c1_1","c2_1","c3_1","c4_1","c5_1", "c1_2","c2_2","c3_2","c4_2","c5_2")
  WMExtExpected <- as(WMExtExpected, "dgCMatrix")
  
  saveit( WMExtExpected=WMExtExpected,
          file = paste0(network_dir, "WMExtExpected.rda"))
}
test.equality <- function(WM_ID, isExtended=FALSE){

  if (!isExtended){
    load(paste0(network_dir, "WMExpected.rda"))
  }else{
    load(paste0(network_dir, "WMExtExpected.rda"))
  } 
  
  load(paste0(network_dir, "WM_", WM_ID, ".rda"))

  print("Col sums for WM")
  print(colSums(WM))  
  test_that("Check that Matrices computation is ok", {
    if(isExtended){
      expect_equal(WM, WMExtExpected, tolerance = 0.00001)
    }else{
      expect_equal(WM, WMExpected, tolerance = 0.00001)
    }
  })
  print("***************         Comparison is SUCCESSFUL         ***************")
}
test.wm.colsums <- function(network_dir, WM_ID){
  load(paste0(network_dir, "WM_", WM_ID, ".rda"))
  test_that("Check that Matrix col sums is ok", {
    expect_equal(unname(colSums(WM)), rep(1, ncol(WM)), tolerance = 0.0001)
  })
  print("Successful: All col sums of WM are 1!")
  WM
}
compare.multiple.matrices <- function(network_dir, WMIDs){
  for(i in WMIDs){
    isEqual <- TRUE
    print(paste0("-------> Testing for WM - ", i))
    load(paste0(network_dir, "WM_", i, ".rda"))
    WM_new <- WM
    load(paste0(network_dir,"WM_", 1, ".rda"))
    WM_old <- WM
    
    if (dim(WM_new)[1]<=46000){
      if (all(WM_new!=WM_old) ) isEqual <- FALSE
    }else{
      print(paste0("WM is too large for lumpsum comparison! Size is ",dim(WM_new)[1]))
      st <- seq(1, dim(WM_new)[1], 46000)
      en <- seq(46000, dim(WM_new)[1], 46001)
      if(!dim(WM_new)[1]%in%en) en <- c(en, dim(WM_new)[1])
      for(j in 1:length(st)){
        print(paste0("Testing for indices : ", st[j], " - ", en[j]))
        if (all(WM_new[st[j]:en[j],st[j]:en[j]]!=WM_old[st[j]:en[j],st[j]:en[j]])){
          isEqual <- FALSE
          break
        }
      }
    }
    if(isEqual)  print(paste0("-------> WM - ", i, " successful"))
    else  print(paste0("-------> WM - ", i, " failed!"))
    cat("\n")
  }
}

create.extended.test.networks()
# create.base.test.networks()
test.equality(1, isExtended = T)

#compare.multiple.matrices(network_dir, c(5, 6))

