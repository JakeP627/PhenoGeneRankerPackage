% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PhenoGeneRankerFunctions.R
\name{CreateWalkMatrix}
\alias{CreateWalkMatrix}
\title{Create Walk Matrix}
\usage{
CreateWalkMatrix(
  inputFileName,
  numCores = 1,
  delta = 0.5,
  zeta = 0.5,
  lambda = 0.5
)
}
\arguments{
\item{inputFileName}{The name of the text file that contains the name of gene and phenotype data files. Has to be a '.txt' file. This text file is made from two columns that are tab-separated. the first row needs the two columns to say type and file_name. Every row after that contains what kind of data each file would contain. These can be either gene, phenotype, or bipartite files. Under file_name, you will input the  name of the file along with the .txt extension. Inside each file_name file, there is a from, to, and weight column which are all tab-separated. For gene and phenotype layer files, the from and to columns will have genes and phenotypes, the order doesn’t matter. The weight column will have the weight of the relationship of the genes and phenotypes, for unweighted network all weights must have value of 1. For a bipartite file, the from column must have genes, the to column must have phenotypes, the meaning and usage of weight is similar to the gene and phenotype layers.}

\item{numCores}{This is the number of cores used for parallel processing.}

\item{delta}{This is the probability of jumping between gene layers. It has a range of 0 to 1. It has a default value of 0.5.}

\item{zeta}{This is the probability of jumping between gene layers. It has a range of 0 to 1. It has a default value of 0.5.}

\item{lambda}{This is the Inter-Network Jump Probability. It has a range of 0 to 1. It has a default value of 0.5.}
}
\value{
This returns a list that contains a walk matrix, list of sorted genes by name in the multiplex heterogeneous network, list of phenotypes sorted by name in the multiplex heterogeneous networkpool nodes, the degree of all genes in the multiplex heterogeneous network to be used for p-value generation in randomWalkRestarts method, the degree of all phenotypes in the multiplex heterogeneous network to be used for p-value generation in randomWalkRestarts method, the number of gene layer files there are, the number of phenotype layer files there are, the number of gene pool nodes, and the number of phenotype nodes.
}
\description{
Generates a Walk matrix from Gene and Phenotype data.
}
\examples{
CreateWalkMatrix('myInput.txt')
CreateWalkMatrix('file.txt', detectCores(), 0.4, 0.7, 0.9)
}
