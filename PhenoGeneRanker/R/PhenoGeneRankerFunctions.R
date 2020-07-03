read.settings <- function(input.file){

  #Layers <- vector("list")
  #Layers[["type"]] <- c("ppi", "pwy", "expr","pheno")
  #Layers[["color"]] <- c("blue", "green", "red","darkorchid")
  #Layers[["count"]] <- c(rep(0, length(Layers[["type"]])))

  files <-  read.table(input.file, header=TRUE,sep="\t", stringsAsFactors = FALSE)


  #Parameters_File <- read.csv(files$file_name[files$layer_type=="parameters"], header=TRUE,sep="\t",dec=".",stringsAsFactors = FALSE)

  LG <- length(files$layer_name[files$type=="gene"])
  LP <- length(files$layer_name[files$type=="phenotype"])
  #Parameters <- check.parameters(Parameters_File, LG, LP)

  output=list(
    FilesDF=files,
    LG=LG,
    LP=LP)
}
read.network.layers <- function(filesDF){
  # read the gene layer names
  filesDF <- filesDF[filesDF$type%in%c("gene", "phenotype", "bipartite"),]
  NetworkLayers <- vector("list",nrow(filesDF))
  j <- 1
  for(f in filesDF$file_name){
    NetworkLayers[[j]][["DF"]] <-  read.table(f, header=TRUE, sep="\t", stringsAsFactors = FALSE)
    NetworkLayers[[j]][["type"]] <- filesDF$type[j]
    
    NetworkLayers[[j]][["graph"]]  <- graph.data.frame(NetworkLayers[[j]][["DF"]], directed = FALSE)
    
    NetworkLayers[[j]][["name"]] <- f
    j <- j+1
  }
  
  gene_pool_nodes_sorted <- generate.pool.nodes(NetworkLayers, type = "gene")
  phenotype_pool_nodes_sorted <- generate.pool.nodes(NetworkLayers, type = "phenotype")
  
  
  idx <- which(lapply(NetworkLayers, `[[`, "type") == "gene")
  for(i in idx){
    NetworkLayers[[i]][["graph"]] <- add.missing.nodes.to.graph(NetworkLayers[[i]][["graph"]],
                                                                "gene", gene_pool_nodes_sorted)  
  }
  
  idx <- which(lapply(NetworkLayers, `[[`, "type") == "phenotype")
  for(i in idx){
    NetworkLayers[[i]][["graph"]] <- add.missing.nodes.to.graph(NetworkLayers[[i]][["graph"]],
                                                                "phenotype", phenotype_pool_nodes_sorted)  
  }  
  
  output=list(
    NetworkLayers=NetworkLayers,
    gene_pool_nodes_sorted=gene_pool_nodes_sorted,
    phenotype_pool_nodes_sorted=phenotype_pool_nodes_sorted)    
  return(output)
}

read.seeds <- function(fileName, gene_pool_nodes_sorted, pheno_pool_nodes_sorted){
  AllSeeds <- read.csv(fileName, header=FALSE, sep="\t", stringsAsFactors = FALSE)
  AllSeeds <- AllSeeds$V1
  SeedList <- check.seeds(AllSeeds,gene_pool_nodes_sorted,pheno_pool_nodes_sorted)
  return(SeedList)
}
generate.pool.nodes <- function(FullNet, type){
  idx <- which(lapply(FullNet, `[[`, "type") == type)
  DFs <- lapply(FullNet[idx], `[[`, "DF")
  Node_Names_all <- unique(c(unlist(lapply(DFs, '[[', 'from')), unlist(lapply(DFs, '[[', 'to'))))

  ## We remove duplicates and sort
  pool_nodes_sorted <- sort(Node_Names_all)

  return(pool_nodes_sorted)
}
add.missing.nodes.to.graph <- function(g, type, pool_nodes_sorted){
  ## We add to each layer the missing nodes of the total set of nodes, of the pool of nodes.
  Node_Names_Layer <- V(g)$name
  Missing_Nodes <- pool_nodes_sorted[which(!pool_nodes_sorted %in% Node_Names_Layer)]
  g <- add_vertices(g ,length(Missing_Nodes), name=Missing_Nodes)
  return(g)
}

#' This method allows you to import all the parameters that you need for use with ease of use for the other methods.
#'
#' @param delta This is the probability of jumping between gene layers. It has a range of 0 to 1.
#' @param zeta This is the probability of jumping between gene layers. It has a range of 0 to 1.
#' @param lambda This is the Inter Network Jump Probability. It has a range of 0 to 1.
#' @param weight This is the weight given to the matrix. it has a range from 0 to 1.
#'
#' @return This returns a list of the parameters.
#' @export
#'
#' @examples
addParameters <- function(delta, zeta, lambda, weight){

  if (delta > 1 || delta< 0){ stop("Incorrect delta, it must be between 0 and 1")}

  if (zeta > 1 || zeta < 0){ stop("Incorrect zeta, it must be between 0 and 1")}

  if (lambda > 1 || lambda < 0){ stop("Incorrect lambda, it must be between 0 and 1")}

  if (weight != 1 & weight != 0){ stop("Incorrect weighted, it must be either 0 or 1")}

  parameters <- list(delta, zeta, lambda, weight)
  names(parameters) <- c("delta", "zeta", "lambda", "weighted")
  return(parameters)
}
check.seeds <- function(Seeds, All_genes,All_Phenotypes){

  Genes_Seeds_Ok <- Seeds[which(Seeds %in% All_genes)]
  Phenotype_Seeds_Ok <- Seeds[which(Seeds %in% All_Phenotypes)]
  All_seeds_ok <- c(Genes_Seeds_Ok,Phenotype_Seeds_Ok)
  All_seeds_ko <- Seeds[which(!Seeds %in% All_seeds_ok)]

  list_Seeds_Ok <- list(Genes_Seeds_Ok,Phenotype_Seeds_Ok)

  print("Seeds OK: ")
  print(paste(All_seeds_ok, sep=" "))
  print("Seeds KO: ")
  print(paste(All_seeds_ko, sep=" "))

  if ((length(Genes_Seeds_Ok) == 0) &&  (length(Phenotype_Seeds_Ok) ==0)){
    stop("Seeds not found in our network")
  } else {
    return(list_Seeds_Ok)
  }

}

get.seed.scores <- function(GeneSeeds, PhenoSeeds, eta, LG, LP, tau, phi) {

  n <- length(GeneSeeds)
  m <- length(PhenoSeeds)

  if ((n != 0 && m!= 0)){

    Seed_Genes_Layer_Labeled <- paste0(rep(GeneSeeds,LG), sep="_",rep(seq(LG), length.out = n*LG,each=n))
    Seeds_Genes_Scores <- rep(((1-eta) * tau)/n,n)

    Seed_Phenos_Layer_Labeled <- paste0(rep(PhenoSeeds,LP), sep="_",rep(seq(LP), length.out = m*LP,each=m))
    Seeds_Phenos_Scores <- rep((eta * phi)/m,m)

  } else {
    eta <- 1
    if (n == 0){
      Seed_Genes_Layer_Labeled <- character()
      Seeds_Genes_Scores <- numeric()

      Seed_Phenos_Layer_Labeled <- paste0(rep(PhenoSeeds,LP), sep="_",rep(seq(LP), length.out = m*LP,each=m))
      Seeds_Phenos_Scores <- rep((eta * phi)/m, m)
    } else {
      Seed_Genes_Layer_Labeled <- paste0(rep(GeneSeeds,LG), sep="_",rep(seq(LG), length.out = n*LG,each=n))
      Seeds_Genes_Scores <- rep(tau/n, n)

      Seed_Phenos_Layer_Labeled <- character()
      Seeds_Phenos_Scores <- numeric()
    }
  }

  ### We prepare a data frame with the seeds.
  Seeds_Score <- data.frame(Seeds_ID = c(Seed_Genes_Layer_Labeled, Seed_Phenos_Layer_Labeled),
                            Score = c(Seeds_Genes_Scores, Seeds_Phenos_Scores)  ,stringsAsFactors = FALSE)
  return(Seeds_Score)
}
categorize.ranked.genes <- function(CandidateGenes, Final_Rank_Genes, generate.p_values=TRUE){
  if(generate.p_values){
    Final_Rank_Genes$Rank <- as.numeric(rownames(Final_Rank_Genes))
  }
  Final_Rank_Genes$Type <- "_"
  Final_Rank_Genes$Type[Final_Rank_Genes$Gene %in% CandidateGenes] <- "C"
  # reorder so that 3rd column will be Type
  #Final_Rank_Genes <- Final_Rank_Genes[,c(1, 2, 4, 3 )]
  Final_Rank_Genes <- Final_Rank_Genes[,c(1, (ncol(Final_Rank_Genes)-1), ncol(Final_Rank_Genes), 2: (ncol(Final_Rank_Genes)-2))]
  Final_Rank_Genes
}
create.supraadjacency.matrix <- function(WholeNet, type, N, L, zeta, is.weighted.graph){
  # WholeNet <- FullNet
  # type="gene"
  # L <- LG
  # zeta <- Parameters$zeta
  # is.weighted.graph <- TRUE
  Graphs <- WholeNet[which(lapply(WholeNet, `[[`, "type") == type)]
  Graphs <- lapply(Graphs, `[[`, "graph")

  Idem_Matrix <- Diagonal(N, x = 1)

  Col_Node_Names <- character()
  Row_Node_Names <- character()

  t1 <- Sys.time()
  if (getDoParWorkers()>1) registerDoSEQ()
  cl <- makeCluster(L+1)
  registerDoParallel(cl)
  SupraAdjacencyResult <- foreach (i=1:L, .packages=c('igraph', 'Matrix')) %dopar% {

    SupraAdjacencyMatrixPart <- Matrix(0,ncol=N*L,nrow=N,sparse = TRUE)

    Adjacency_Layer <-  as_adjacency_matrix(Graphs[[i]], attr="weight", sparse = TRUE)

    if (!is.weighted.graph){
      Adjacency_Layer <-  as_adjacency_matrix(Graphs[[i]], sparse = TRUE)
    }

    ## We order the matrix by the node name. This way all the matrix will have the same. Additionally we include a label with the layer number for each node name.
    Adjacency_Layer <- Adjacency_Layer[order(rownames(Adjacency_Layer)),order(colnames(Adjacency_Layer))]
    Layer_Row_Col_Names <- paste(colnames(Adjacency_Layer),i,sep="_")

    ## We fill the diagonal blocks with the adjacencies matrix of each layer.
    Position_ini_row <- 1
    Position_end_row <- N
    Position_ini_col <- 1 + (i-1)*N
    Position_end_col <- N + (i-1)*N
    SupraAdjacencyMatrixPart[(Position_ini_row:Position_end_row),(Position_ini_col:Position_end_col)] <- (1-zeta)*(Adjacency_Layer)

    # avoid division by zero for monoplex network
    L_mdfd <- L-1
    if (L == 1) L_mdfd <- 1

    ## We fill the off-diagonal blocks with the transition probability among layers.
    for (j in 1:L){
      Position_ini_col <- 1 + (j-1)*N
      Position_end_col <- N + (j-1)*N
      if (j != i){
        SupraAdjacencyMatrixPart[(Position_ini_row:Position_end_row),(Position_ini_col:Position_end_col)] <- (zeta/(L_mdfd))*Idem_Matrix
      }
    }
    return(list(SupraAdjacencyMatrixPart, Layer_Row_Col_Names))
  }

  stopCluster(cl)
  #cat(type, " SupraAdjacency Time for parallel part: ", format(Sys.time()-t1), "\n")
  t2 <- Sys.time()
  #SupraAdjacencyMatrix <- Matrix(0,ncol=N*L,nrow=N*L,sparse = TRUE)
  SupraAdjacencyResult <- unlist(SupraAdjacencyResult, recursive = FALSE)


  # Row-Col names are even indexed
  Col_Names <- do.call('c',SupraAdjacencyResult[seq(2,2*L,by=2)])

  # Parallele parts of the SupraAdjacencyMatrix are odd indexed
  SupraAdjacencyMatrix <- do.call('rbind', SupraAdjacencyResult[seq(1,2*L,by=2)])

  #SupraAdjacencyMatrix <- rbind.fill.matrix(SupraAdjacencyResult[seq(1,2*L,by=2)])

  #cat("SupraAdjacency Time for rest: ", format(Sys.time()-t2), "\n")
  rownames(SupraAdjacencyMatrix) <- Col_Names
  colnames(SupraAdjacencyMatrix) <- Col_Names

  return(SupraAdjacencyMatrix)
}
create.bipartite.matrix <- function(WholeNet, N, M, gene_pool_nodes_sorted, pheno_pool_nodes_sorted, numCores, isWeighted){
  #WholeNet <- FullNet
  Gene_Phenoivar_Network <- WholeNet[which(lapply(WholeNet, `[[`, "type") == "bipartite")]
  Gene_Phenoivar_Network <- lapply(Gene_Phenoivar_Network, `[[`, "DF")[[1]]

  # Get the Subset of Gene-Phenoivar relations which have common genes in whole network
  Gene_Phenoivar_Network <- Gene_Phenoivar_Network[which(Gene_Phenoivar_Network$from %in% gene_pool_nodes_sorted), ]
  Gene_Phenoivar_Network <- Gene_Phenoivar_Network[which(Gene_Phenoivar_Network$to %in% pheno_pool_nodes_sorted), ]

  Gene_Phenoivar_Network <- graph.data.frame(Gene_Phenoivar_Network, directed = FALSE)

  el <- as_edgelist(Gene_Phenoivar_Network)
  value <- edge_attr(Gene_Phenoivar_Network, name = "weight")
  if (!isWeighted){
    value <- rep(1, nrow(el))
  }

  Bipartite_matrix <- Matrix(data=0, nrow=N, ncol=M)
  rownames(Bipartite_matrix) <- gene_pool_nodes_sorted
  colnames(Bipartite_matrix) <- pheno_pool_nodes_sorted
  # rindx <- unlist(mclapply(el[,1], function(x) which(rownames(Bipartite_matrix) %in% x), mc.cores=10))
  # cindx <- unlist(mclapply(el[,2], function(x) which(colnames(Bipartite_matrix) %in% x), mc.cores=10))
  rindx <- unlist(lapply(el[,1], function(x) which(rownames(Bipartite_matrix) %in% x)))
  cindx <- unlist(lapply(el[,2], function(x) which(colnames(Bipartite_matrix) %in% x)))

  lenRindx <- length(rindx)
  partLen <- floor(lenRindx/numCores)

  if (partLen == 0){
    cat("numCores:", numCores, " - partLen: ", partLen, "\n")
    stop("Assigned numCores is greater than data length! Assign less!")
  }

  rindx_parts <- list()
  cindx_parts <- list()
  value_parts <- list()
  for(i in 1:numCores){
    stInd <- (i-1)*partLen + 1
    endInd <- i*partLen
    if (i==numCores){
      endInd <- lenRindx
    }
    rindx_parts[[i]] <- rindx[stInd:endInd]
    cindx_parts[[i]] <- cindx[stInd:endInd]
    value_parts[[i]] <- value[stInd:endInd]
  }


  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  Bipartite_matrix_result <- foreach (i=1:numCores, .packages='Matrix') %dopar% {
    for(j in 1:length(rindx_parts[[i]])){
      Bipartite_matrix[rindx_parts[[i]][j],cindx_parts[[i]][j]] <- value_parts[[i]][j]
    }
    return(Bipartite_matrix)
  }
  stopCluster(cl)
  Bipartite_matrix <- Reduce('+', Bipartite_matrix_result)

  return(Bipartite_matrix)
}
create.suprabipartite.matrix <- function(Bipartite_matrix, N, M, LG, LP){
  SupraBipartiteMatrix <- Matrix(0,nrow=N*LG,ncol=M*LP,sparse = TRUE)

  Row_Node_Names <- sprintf(paste0(rep(rownames(Bipartite_matrix), LG), "_%d"),
                            rep(seq_len(LG), each=N))
  SupraBipartiteMatrix <- do.call(rbind, replicate(LG, Bipartite_matrix, simplify=FALSE))

  rownames(SupraBipartiteMatrix) <- Row_Node_Names


  Col_Node_Names <- sprintf(paste0(rep(colnames(Bipartite_matrix), LP), "_%d"),
                            rep(seq_len(LP), each=M))

  SupraBipartiteMatrix <- do.call(cbind, replicate(LP, SupraBipartiteMatrix, simplify=FALSE))
  colnames(SupraBipartiteMatrix) <- Col_Node_Names
  return(SupraBipartiteMatrix)
}
create.transition.matrix <- function(SupraBipartiteMatrix, N, M, LG, LP, lambda, isTranspose){
  ###TRUE = Row wise/pheno-gene, FALSE = Col Wise/gene-pheno
  if(isTranspose){
    #### Transition Matrix for the inter-subnetworks links
    TransitionMatrix <- Matrix(0,nrow=M*LP,ncol=N*LG,sparse = TRUE)
    colnames(TransitionMatrix) <- rownames(SupraBipartiteMatrix)
    rownames(TransitionMatrix) <- colnames(SupraBipartiteMatrix)

    Row_Sum_Bipartite <- rowSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,sparseResult = FALSE)
    # row wise normalization for propability
    for (i in 1:(N*LG)){
      if (Row_Sum_Bipartite[i] != 0){
        TransitionMatrix[,i] <- (lambda*SupraBipartiteMatrix[i,])/Row_Sum_Bipartite[i]
      }
    }
  }
  else{
    #### Transition Matrix for the inter-subnetworks links
    TransitionMatrix <- Matrix(0,nrow=N*LG,ncol=M*LP,sparse = TRUE)
    colnames(TransitionMatrix) <- colnames(SupraBipartiteMatrix)
    rownames(TransitionMatrix) <- rownames(SupraBipartiteMatrix)

    Col_Sum_Bipartite <- colSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,sparseResult = FALSE)

    # columnwise normalization for propability
    for (j in 1:(M*LP)){
      if (Col_Sum_Bipartite[j] != 0){
        TransitionMatrix[,j] <- (lambda*SupraBipartiteMatrix[,j]) /Col_Sum_Bipartite[j]
      }
    }
  }


  return(TransitionMatrix)
}
create.transition.matrix.gene_pheno <- function(SupraBipartiteMatrix, N, M, LG, LP, lambda){
  #### Transition Matrix for the inter-subnetworks links
  Transition_Gene_Phenoivar <- Matrix(0,nrow=N*LG,ncol=M*LP,sparse = TRUE)
  colnames(Transition_Gene_Phenoivar) <- colnames(SupraBipartiteMatrix)
  rownames(Transition_Gene_Phenoivar) <- rownames(SupraBipartiteMatrix)

  Col_Sum_Bipartite <- colSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,sparseResult = FALSE)

  # columnwise normalization for propability
  for (j in 1:(M*LP)){
    if (Col_Sum_Bipartite[j] != 0){
      Transition_Gene_Phenoivar[,j] <- (lambda*SupraBipartiteMatrix[,j]) /Col_Sum_Bipartite[j]
    }
  }
  return(Transition_Gene_Phenoivar)
}
create.transition.matrix.pheno_gene <- function(SupraBipartiteMatrix, N, M, LG, LP, lambda){
  Transition_Phenoivar_Gene <- Matrix(0,nrow=M*LP,ncol=N*LG,sparse = TRUE)

  colnames(Transition_Phenoivar_Gene) <- rownames(SupraBipartiteMatrix)
  rownames(Transition_Phenoivar_Gene) <- colnames(SupraBipartiteMatrix)

  Row_Sum_Bipartite <- rowSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,sparseResult = FALSE)
  # row wise normalization for propability
  for (i in 1:(N*LG)){
    if (Row_Sum_Bipartite[i] != 0){
      Transition_Phenoivar_Gene[,i] <- (lambda*SupraBipartiteMatrix[i,])/Row_Sum_Bipartite[i]
    }
  }
  return(Transition_Phenoivar_Gene)
}
create.transition.multiplex.network <- function(SupraAdjacencyMatrix, SupraBipartiteMatrix, Num, inputLength, lambda, numCores){
  #### Transition Matrix for the intra-subnetworks links
  Transition_Multiplex_Network <- Matrix(0,nrow=Num*inputLength,ncol=Num*inputLength,sparse = TRUE)

  rownames(Transition_Multiplex_Network) <- rownames(SupraAdjacencyMatrix)
  colnames(Transition_Multiplex_Network) <- colnames(SupraAdjacencyMatrix)

  Col_Sum_Multiplex <- colSums(SupraAdjacencyMatrix,na.rm=FALSE,dims=1, sparseResult=FALSE)
  Row_Sum_Bipartite <- rowSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,sparseResult = FALSE)


  partLen <- floor(Num*inputLength/numCores)
  # below can only happen with toy samples
  if (partLen==0){
    stop("Assigned numCores is greater than data length! Assign less!")
  }

  parts <- vector('list',numCores)
  for(i in 1:numCores){
    stInd <- (i-1)*partLen + 1
    endInd <- i*partLen
    if (i==numCores){
      endInd <- Num*inputLength
    }
    parts[[i]][['start']] <- stInd
    parts[[i]][['end']] <- endInd
  }

  cl <- makeCluster(numCores+1)

  registerDoParallel(cl)

  Transition_Multiplex_Network_Result <- foreach (i=1:numCores, .packages='Matrix') %dopar% {
    for (j in parts[[i]]["start"]:parts[[i]]["end"]){
      if(Row_Sum_Bipartite[j] != 0){
        Transition_Multiplex_Network[,j] <- ((1-lambda)*SupraAdjacencyMatrix[,j]) /Col_Sum_Multiplex[j]
      } else {
        Transition_Multiplex_Network[,j] <- SupraAdjacencyMatrix[,j] /Col_Sum_Multiplex[j]
      }
    }
    return(Transition_Multiplex_Network)
  }

  Transition_Multiplex_Network <- Reduce('+', Transition_Multiplex_Network_Result)
  stopCluster(cl)

  return(Transition_Multiplex_Network)
}
create.gene.transition.multiplex.network <- function(SupraAdjacencyMatrix, SupraBipartiteMatrix, N, LG, lambda, numCores){
  #### Transition Matrix for the intra-subnetworks links
  Transition_Multiplex_Network <- Matrix(0,nrow=N*LG,ncol=N*LG,sparse = TRUE)

  rownames(Transition_Multiplex_Network) <- rownames(SupraAdjacencyMatrix)
  colnames(Transition_Multiplex_Network) <- colnames(SupraAdjacencyMatrix)

  Col_Sum_Multiplex <- colSums(SupraAdjacencyMatrix,na.rm=FALSE,dims=1, sparseResult=FALSE)
  Row_Sum_Bipartite <- rowSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,sparseResult = FALSE)


  partLen <- floor(N*LG/numCores)
  # below can only happen with toy samples
  if (partLen==0){
    stop("Assigned numCores is greater than data length! Assign less!")
  }

  parts <- vector('list',numCores)
  for(i in 1:numCores){
    stInd <- (i-1)*partLen + 1
    endInd <- i*partLen
    if (i==numCores){
      endInd <- N*LG
    }
    parts[[i]][['start']] <- stInd
    parts[[i]][['end']] <- endInd
  }

  cl <- makeCluster(numCores+1)

  registerDoParallel(cl)

  Transition_Multiplex_Network_Result <- foreach (i=1:numCores, .packages='Matrix') %dopar% {
    for (j in parts[[i]]["start"]:parts[[i]]["end"]){
      if(Row_Sum_Bipartite[j] != 0){
        Transition_Multiplex_Network[,j] <- ((1-lambda)*SupraAdjacencyMatrix[,j]) /Col_Sum_Multiplex[j]
      } else {
        Transition_Multiplex_Network[,j] <- SupraAdjacencyMatrix[,j] /Col_Sum_Multiplex[j]
      }
    }
    return(Transition_Multiplex_Network)
  }

  Transition_Multiplex_Network <- Reduce('+', Transition_Multiplex_Network_Result)
  stopCluster(cl)

  return(Transition_Multiplex_Network)
}
create.pheno.transition.multiplex.network <- function(SupraAdjacencyMatrixPheno, SupraBipartiteMatrix, M, LP, lambda, numCores){
  Transition_Multiplex_Network_Pheno <- Matrix(0,nrow=M*LP,ncol=M*LP,sparse = TRUE)

  rownames(Transition_Multiplex_Network_Pheno) <- rownames(SupraAdjacencyMatrixPheno)
  colnames(Transition_Multiplex_Network_Pheno) <- colnames(SupraAdjacencyMatrixPheno)

  Col_Sum_Multiplex <- colSums(SupraAdjacencyMatrixPheno,na.rm=FALSE,dims=1, sparseResult=FALSE)
  Col_Sum_Bipartite <- colSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,sparseResult = FALSE)

  partLen <- floor(M*LP/numCores)
  # below can only happen with toy samples
  if (partLen==0){
    stop("Assigned numCores is greater than data length! Assign less!")
  }

  parts <- vector('list',numCores)
  for(i in 1:numCores){
    stInd <- (i-1)*partLen + 1
    endInd <- i*partLen
    if (i==numCores){
      endInd <- M*LP
    }
    parts[[i]][['start']] <- stInd
    parts[[i]][['end']] <- endInd
  }
  cl <- makeCluster(numCores+1)

  registerDoParallel(cl)
  Transition_Multiplex_Network_Pheno_Result <- foreach (i=1:numCores, .packages='Matrix') %dopar% {
    for (j in parts[[i]]["start"]:parts[[i]]["end"]){
      if(Col_Sum_Bipartite[j] != 0){
        Transition_Multiplex_Network_Pheno[,j] <- ((1-lambda)*SupraAdjacencyMatrixPheno[,j]) /Col_Sum_Multiplex[j]
      } else {
        Transition_Multiplex_Network_Pheno[,j] <- SupraAdjacencyMatrixPheno[,j] /Col_Sum_Multiplex[j]
      }
    }
    return(Transition_Multiplex_Network_Pheno)
  }

  Transition_Multiplex_Network_Pheno <- Reduce('+', Transition_Multiplex_Network_Pheno_Result)
  stopCluster(cl)
  return(Transition_Multiplex_Network_Pheno)
}

rank_genes <- function(Number_Genes, Number_Layers,Results,Seeds){
  ## We sort the score to obtain the ranking of Genes and Phenotypes.
  genes_rank <- data.frame(Gene = character(length = Number_Genes), Score = 0)
  genes_rank$Gene <- gsub("_1", "", row.names(Results)[1:Number_Genes])

  ## We calculate the Geometric Mean among the genes in the different layers.
  genes_rank$Score <- geometric.mean(as.vector(Results[,1]),Number_Layers,Number_Genes)

  genes_rank_sort <- genes_rank[with(genes_rank, order(-Score, Gene)), ]

  ### We remove the seed genes from the Ranking
  genes_rank_sort_NoSeeds <- genes_rank_sort[which(!genes_rank_sort$Gene %in% Seeds),]

  genes_rank_sort_NoSeeds$Rank <- seq(1, nrow(genes_rank_sort_NoSeeds))
  genes_rank_sort_NoSeeds <-genes_rank_sort_NoSeeds[,c("Rank","Gene","Score")]

  return(genes_rank_sort_NoSeeds)
}
rank_phenotypes <- function(Number_Genes,Num_Gene_Layers,Number_Phenotypes,Num_Pheno_Layers, Results,Seeds){

  ## rank_phenotypes
  phenotypes_rank <- data.frame(Pheno = character(length = Number_Phenotypes), Score = 0)
  phenotypes_rank$Pheno <- gsub("_1", "", row.names(Results)[(Number_Genes*Num_Gene_Layers+1):(Number_Genes*Num_Gene_Layers+Number_Phenotypes)])

  phenotypes_rank$Score <- geometric.mean(as.vector(Results[,1])[(Number_Genes*Num_Gene_Layers+1):nrow(Results)],Num_Pheno_Layers,Number_Phenotypes)

  phenotypes_rank_sort <- phenotypes_rank[with(phenotypes_rank, order(-Score, Pheno)), ]
  phenotypes_rank_sort_NoSeeds <- phenotypes_rank_sort[which(!phenotypes_rank_sort$Pheno %in% Seeds),]

  phenotypes_rank_sort_NoSeeds$Rank <- seq(1, nrow(phenotypes_rank_sort_NoSeeds))
  phenotypes_rank_sort_NoSeeds <-phenotypes_rank_sort_NoSeeds[,c("Rank","Pheno","Score")]

  return(phenotypes_rank_sort_NoSeeds)
}
geometric.mean <- function(Scores, L, N) {

  FinalScore <- numeric(length = N)

  for (i in seq_len(N)){
    FinalScore[i] <- prod(Scores[seq(from = i, to = N*L, by=N)])^(1/L)
  }

  return(FinalScore)
}

#' Title
#'
#' @param Walk_Matrix 
#' @param GeneSeeds 
#' @param PhenoSeeds 
#' @param eta 
#' @param LG 
#' @param LP 
#' @param tau 
#' @param phi 
#'
#' @return
#' @export
#'
#' @examples
randomWalkRestarts <- function(Walk_Matrix, GeneSeeds, PhenoSeeds, LG, LP, tau=FALSE, phi=FALSE, eta=0.5){
  ### We define the threshold and the number maximum of iterations for the randon walker.
  if(tau==FALSE){
    tau <- rep(1,LG)
  }
  if(phi==FALSE){
    phi <- rep(1,LG)
  }
    
  Seeds_Score <- get.seed.scores(GeneSeeds,PhenoSeeds, eta, LG, LP, Parameters$tau/LG, Parameters$phi/LP)
  Threeshold <- 1e-10
  NetworkSize <- ncol(Walk_Matrix)

  ### We initialize the variables to control the flux in the RW algo.
  residue <- 1
  iter <- 1

  #### We define the prox_vector(The vector we will move after the first RW iteration. We start from The seed. We have to take in account
  #### that the walker with restart in some of the Seed genes, depending on the score we gave in that file).
  prox_vector <- Matrix(0,nrow = NetworkSize,ncol=1, sparse=TRUE)

  prox_vector[which(colnames(Walk_Matrix) %in% Seeds_Score[,1])] <- (Seeds_Score[,2])

  prox_vector  <- prox_vector/sum(prox_vector)
  restart_vector <-  prox_vector
  while(residue >= Threeshold){
    old_prox_vector <- prox_vector
    prox_vector <- (1-r)*(Walk_Matrix %*% prox_vector) + r*restart_vector

    residue <- sqrt(sum((prox_vector-old_prox_vector)^2))

    iter <- iter + 1;
  }

   print("RWR-MH number of iteration: ")
   print(iter-1)
  return(prox_vector)
}
get.connectivity <- function(NetworkDF, gene_pool_nodes_sorted, pheno_pool_nodes_sorted){
  WholeNet <- bind_rows(lapply(NetworkDF, `[[`, "DF"))
  g <- graph.data.frame(WholeNet, directed=FALSE)
  A <- as_adjacency_matrix(g, sparse = TRUE, attr="weight")
  Degree=apply(A, 2, sum)
  Connectivity <- data.frame(Node = as.character(A@Dimnames[[1]]), Degree = Degree, row.names = NULL)
  Connectivity <- Connectivity[order(Connectivity$Degree, decreasing = TRUE),]
  Connectivity$Node <- as.character(Connectivity$Node)
  GeneConnectivity <- Connectivity[which(Connectivity$Node %in% gene_pool_nodes_sorted),]
  PhenoConnectivity <- Connectivity[which(Connectivity$Node %in% pheno_pool_nodes_sorted),]
  return(list(gene=GeneConnectivity,pheno=PhenoConnectivity))
}


get.WalkMatrixID <- function(filesDF){
  WM_DF <-  read.table(filesDF$file_name[filesDF$layer_type=="WM_DF"], header=TRUE, sep="\t", stringsAsFactors = FALSE)
  LayerNames <- filesDF$layer_name[filesDF$type%in%c("gene", "pheno", "bipartite")]
  LayerNames <- paste0(LayerNames,collapse=',')
  new_WM <- data.frame(ID=ifelse(length(WM_DF$ID)==0, 1, max(WM_DF$ID)+1),Layers=LayerNames, Date=as.character(Sys.time()))
  WM_DF <- rbind(WM_DF, new_WM)
  res <- try(write.table(WM_DF, filesDF$file_name[filesDF$layer_type=="WM_DF"], sep="\t",
              col.names = TRUE, row.names = FALSE, quote=FALSE))
  if(!is.null(res)) stop("Eror: Could not get Walk Matrix ID!")
  output=list(
    ID=max(WM_DF$ID),
    LayerNames=LayerNames)
  return(output)
}


#' Generates a Walk matrix file from Gene and Phenoivar data
#'
#' Generates a walk matrix file from the data given by the user.
#' The data for each gene and phenotype is added to supraadjency matrixes,
#' bipartite matrixes, and suprabipartite matrixes.These are combined to make transition and
#' transistion multiplex networks. A multiplex heterogeneous network is created by
#' combining the transition and transition multiplex networks.
#' That network then gets saved to a list
#'
#' @param inputFileName name of the txt file that contains the gene and phenotype data.
#' @param numCores number of cores used for parallel processing.
#'     - recommended to use detectCores() method from the parallel package.
#'
#' @return This returns a walk matrix as a list.
#' @export
#'
#' @examples
createWalkMatrix <- function(inputFileName, numCores, delta=0.5, zeta=0.5, lambda=0.5, weight=1){
  network_range <<- c(0.001, 1)
  global_t1 <- Sys.time()
  t1 <- Sys.time()
  Settings <- read.settings(inputFileName)
  Parameters <- addParameters(delta, zeta, lambda, weight)
  FilesDF <- Settings$FilesDF
  LG <- Settings$LG
  LP <- Settings$LP

  cat("Creating Walk Matrix For: ", paste0(FilesDF$layer_name[FilesDF$type%in%c("gene", "pheno", "bipartite")], collapse = ","), "\n")
  FullNet <- read.network.layers(FilesDF)
  Layers <- FullNet$Layers
  gene_pool_nodes_sorted <- FullNet$gene_pool_nodes_sorted
  pheno_pool_nodes_sorted <- FullNet$pheno_pool_nodes_sorted
  FullNet <- FullNet$NetworkLayers

  N=length(gene_pool_nodes_sorted)
  M <- length(pheno_pool_nodes_sorted)
  cat("Time to Initialize : ", format(Sys.time()-t1), "\n")


  # CREATE SUPRAADJACENCYMATRIX FOR GENES
  t1 <- Sys.time()
  SupraAdjacencyMatrix <- create.supraadjacency.matrix(FullNet, "gene", N, LG, Parameters$zeta, TRUE)
  cat("Time to create SupraAdjacencyMatrix for Genes : ", format(Sys.time()-t1), "\n")


  # CREATE SUPRAADJACENCYMATRIX FOR CULTIVARS
  t1 <- Sys.time()
  SupraAdjacencyMatrixPheno <- create.supraadjacency.matrix(FullNet, "pheno", M, LP, Parameters$delta, TRUE)
  cat("Time to create SupraAdjacencyMatrix for Phenoivars : ", format(Sys.time()-t1), "\n")

  # CREATE THE BIPARTITE GRAPH
  t1 <- Sys.time()
  BipartiteMatrix <- create.bipartite.matrix(FullNet, N, M, gene_pool_nodes_sorted, pheno_pool_nodes_sorted, numCores, Parameters$weighted)
  cat("Time to create BipartiteMatrix : ", format(Sys.time()-t1), "\n")


  # CREATE SUPRABIPARTITEGRAPH

  ## We expand the biparite graph to fit the multiplex dimensions.
  ## The biparti matrix has now NL x MK
  ## The genes in all the layers have to point to the phenotypes in all layers

  t1 <- Sys.time()
  SupraBipartiteMatrix <- create.suprabipartite.matrix(BipartiteMatrix, N, M, LG, LP)
  cat("Time to create SupraBipartiteMatrix : ", format(Sys.time()-t1), "\n")


  t1 <- Sys.time()
  #Transition_Gene_Phenoivar <- create.transition.matrix.gene_pheno(SupraBipartiteMatrix, N, M, LG, LP, Parameters$lambda)
  Transition_Gene_Phenoivar <- create.transition.matrix(SupraBipartiteMatrix, N, M, LG, LP, Parameters$lambda, FALSE)
  cat("Time to create Transition Matrix for Genes-Phenoivars : ", format(Sys.time()-t1), "\n")

  t1 <- Sys.time()
  #Transition_Phenoivar_Gene <- create.transition.matrix.pheno_gene(SupraBipartiteMatrix, N, M, LG, LP, Parameters$lambda)
  Transition_Phenoivar_Gene <- create.transition.matrix(SupraBipartiteMatrix, N, M, LG, LP, Parameters$lambda, TRUE)
  cat("Time to create Transition Matrix for Phenoivars-Genes : ", format(Sys.time()-t1), "\n")


  t1 <- Sys.time()
  Gene_Transition_Multiplex_Network <- create.gene.transition.multiplex.network(SupraAdjacencyMatrix,
                                                                                SupraBipartiteMatrix,
                                                                                N, LG, Parameters$lambda, numCores)
  #Gene_Transition_Multiplex_Network <- create.transition.multiplex.network(SupraAdjacencyMatrix,
   #                                                                             SupraBipartiteMatrix,
    #                                                                            N, LG, Parameters$lambda, numCores)
  cat("Time to create Gene Transition Multiplex Network : ", format(Sys.time()-t1), "\n")

  #printSpMatrix2(Gene_Transition_Multiplex_Network, col.names = TRUE)


  # CREATE TRANSITION MULTIPLEX NETWORK FOR CULTIVARS
  t1 <- Sys.time()
  Pheno_Transition_Multiplex_Network <- create.pheno.transition.multiplex.network(SupraAdjacencyMatrixPheno,
                                                                                SupraBipartiteMatrix,
                                                                                M, LP, Parameters$lambda, numCores)
  #Pheno_Transition_Multiplex_Network <- create.transition.multiplex.network(SupraAdjacencyMatrixPheno,
   #                                                                            SupraBipartiteMatrix,
    #                                                                           M, LP, Parameters$lambda, numCores)
  cat("Time to create Pheno Transition Multiplex Network : ", format(Sys.time()-t1), "\n")


  #printSpMatrix2(round(Pheno_Transition_Multiplex_Network, 3), col.names = TRUE)
  #colSums(Pheno_Transition_Multiplex_Network)

  # CREATE FINAL TRANSITION MULTIPLEX HETERO MATRIX

  t1 <- Sys.time()

  ### We generate the global transiction matrix and we return it.
  Multiplex_Heterogeneous_Matrix <- rbind(cbind(Gene_Transition_Multiplex_Network, Transition_Gene_Phenoivar),
                                          cbind(Transition_Phenoivar_Gene, Pheno_Transition_Multiplex_Network))

  cat("Time to create Final Multiplex Heterogeneous Network : ", format(Sys.time()-t1), "\n")


  #Extract candidate genes for further Random Walks on this WM
  GenePhenoDF <- lapply(FullNet[which(lapply(FullNet, `[[`, "type") == "bipartite")], `[[`, "DF")[[1]]
  CandidateGenes <- unique(GenePhenoDF$from)
  Connectivity <- get.connectivity(FullNet, gene_pool_nodes_sorted, pheno_pool_nodes_sorted)
  #WM_ID <- get.WalkMatrixID(FilesDF)
  WM <- list(WM = Multiplex_Heterogeneous_Matrix,
                  genes = gene_pool_nodes_sorted,
                  phenotypes = phenotype_pool_nodes_sorted,
                  gene_connectivity = Connectivity[["gene"]],
                  phenotype_connectivity = Connectivity[["phenotype"]],
                  LG = LG,
                  LP = LP,
                  N = N,
                  M = M)
  # saveit(WM=Multiplex_Heterogeneous_Matrix,
  #        GeneConnectivity=Connectivity[["gene"]], PhenoConnectivity=Connectivity[["pheno"]],
  #        LG=LG, LP=LP, M=M, N=N, Parameters=Parameters,
  #        gene_pool_nodes_sorted=gene_pool_nodes_sorted,
  #        pheno_pool_nodes_sorted=pheno_pool_nodes_sorted,
  #        CandidateGenes=CandidateGenes,
  #        WM_ID = WM_ID[["ID"]], WM_Layers=WM_ID[["LayerNames"]],
  #        file = paste0("WM_", WM_ID[["ID"]], ".rda"))

  cat("Time to Start-Finish : ", format(Sys.time()-global_t1), "\n")
  # if there are still some paralllel workers stop and force sequential
  registerDoSEQ()
  return(WM)
}

saveit <- function(..., file) {
  x <- list(...)
  save(list=names(x), file=file, envir=list2env(x))
}

createWalkMatrixUpdated <- function(inputFileName, params, numCores){
  network_range <<- c(0.001, 1)
  global_t1 <- Sys.time()
  t1 <- Sys.time()
  Settings <- read.settings(inputFileName)
  Parameters <- params
  FilesDF <- Settings$FilesDF
  LG <- Settings$LG
  LP <- Settings$LP


  cat("Creating Walk Matrix For: ", paste0(FilesDF$layer_name[FilesDF$type%in%c("gene", "pheno", "bipartite")], collapse = ","), "\n")
  FullNet <- read.network.layers(FilesDF, Parameters, Settings$Layers)
  Layers <- FullNet$Layers
  gene_pool_nodes_sorted <- FullNet$gene_pool_nodes_sorted
  pheno_pool_nodes_sorted <- FullNet$pheno_pool_nodes_sorted
  FullNet <- FullNet$NetworkLayers

  N=length(gene_pool_nodes_sorted)
  M <- length(pheno_pool_nodes_sorted)
  cat("Time to Initialize : ", format(Sys.time()-t1), "\n")


  # CREATE SUPRAADJACENCYMATRIX FOR GENES
  t1 <- Sys.time()
  SupraAdjacencyMatrix <- create.supraadjacency.matrix(FullNet, "gene", N, LG, Parameters$zeta, TRUE)
  cat("Time to create SupraAdjacencyMatrix for Genes : ", format(Sys.time()-t1), "\n")


  # CREATE SUPRAADJACENCYMATRIX FOR CULTIVARS
  t1 <- Sys.time()
  SupraAdjacencyMatrixPheno <- create.supraadjacency.matrix(FullNet, "pheno", M, LP, Parameters$delta, TRUE)
  cat("Time to create SupraAdjacencyMatrix for Phenoivars : ", format(Sys.time()-t1), "\n")

  # CREATE THE BIPARTITE GRAPH
  t1 <- Sys.time()
  BipartiteMatrix <- create.bipartite.matrix(FullNet, N, M, gene_pool_nodes_sorted, pheno_pool_nodes_sorted, numCores, Parameters$weighted)
  cat("Time to create BipartiteMatrix : ", format(Sys.time()-t1), "\n")


  # CREATE SUPRABIPARTITEGRAPH

  ## We expand the biparite graph to fit the multiplex dimensions.
  ## The biparti matrix has now NL x MK
  ## The genes in all the layers have to point to the phenotypes in all layers

  t1 <- Sys.time()
  SupraBipartiteMatrix <- create.suprabipartite.matrix(BipartiteMatrix, N, M, LG, LP)
  cat("Time to create SupraBipartiteMatrix : ", format(Sys.time()-t1), "\n")


  t1 <- Sys.time()
  #Transition_Gene_Phenoivar <- create.transition.matrix.gene_pheno(SupraBipartiteMatrix, N, M, LG, LP, Parameters$lambda)
  Transition_Gene_Phenoivar <- create.transition.matrix(SupraBipartiteMatrix, N, M, LG, LP, Parameters$lambda, FALSE)
  cat("Time to create Transition Matrix for Genes-Phenoivars : ", format(Sys.time()-t1), "\n")

  t1 <- Sys.time()
  #Transition_Phenoivar_Gene <- create.transition.matrix.pheno_gene(SupraBipartiteMatrix, N, M, LG, LP, Parameters$lambda)
  Transition_Phenoivar_Gene <- create.transition.matrix(SupraBipartiteMatrix, N, M, LG, LP, Parameters$lambda, TRUE)
  cat("Time to create Transition Matrix for Phenoivars-Genes : ", format(Sys.time()-t1), "\n")


  t1 <- Sys.time()
  Gene_Transition_Multiplex_Network <- create.gene.transition.multiplex.network(SupraAdjacencyMatrix,
                                                                                SupraBipartiteMatrix,
                                                                                N, LG, Parameters$lambda, numCores)
  #Gene_Transition_Multiplex_Network <- create.transition.multiplex.network(SupraAdjacencyMatrix,
  #                                                                             SupraBipartiteMatrix,
  #                                                                            N, LG, Parameters$lambda, numCores)
  cat("Time to create Gene Transition Multiplex Network : ", format(Sys.time()-t1), "\n")

  #printSpMatrix2(Gene_Transition_Multiplex_Network, col.names = TRUE)


  # CREATE TRANSITION MULTIPLEX NETWORK FOR CULTIVARS
  t1 <- Sys.time()
  Pheno_Transition_Multiplex_Network <- create.pheno.transition.multiplex.network(SupraAdjacencyMatrixPheno,
                                                                                SupraBipartiteMatrix,
                                                                                M, LP, Parameters$lambda, numCores)
  #Pheno_Transition_Multiplex_Network <- create.transition.multiplex.network(SupraAdjacencyMatrixPheno,
  #                                                                            SupraBipartiteMatrix,
  #                                                                           M, LP, Parameters$lambda, numCores)
  cat("Time to create Pheno Transition Multiplex Network : ", format(Sys.time()-t1), "\n")


  #printSpMatrix2(round(Pheno_Transition_Multiplex_Network, 3), col.names = TRUE)
  #colSums(Pheno_Transition_Multiplex_Network)

  # CREATE FINAL TRANSITION MULTIPLEX HETERO MATRIX

  t1 <- Sys.time()

  ### We generate the global transiction matrix and we return it.
  Multiplex_Heterogeneous_Matrix <- rbind(cbind(Gene_Transition_Multiplex_Network, Transition_Gene_Phenoivar),
                                          cbind(Transition_Phenoivar_Gene, Pheno_Transition_Multiplex_Network))

  cat("Time to create Final Multiplex Heterogeneous Network : ", format(Sys.time()-t1), "\n")


  #Extract candidate genes for further Random Walks on this WM
  GenePhenoDF <- lapply(FullNet[which(lapply(FullNet, `[[`, "type") == "bipartite")], `[[`, "DF")[[1]]
  CandidateGenes <- unique(GenePhenoDF$from)
  Connectivity <- get.connectivity(FullNet, gene_pool_nodes_sorted, pheno_pool_nodes_sorted)
  #WM_ID <- get.WalkMatrixID(FilesDF)

  WM <- list("WM"=Multiplex_Heterogeneous_Matrix,Connectivity[["gene"]],"GeneConnectivity"=Connectivity[["pheno"]],
             "gene_pool_nodes_sorted"=gene_pool_nodes_sorted, "pheno_pool_nodes_sorted"=pheno_pool_nodes_sorted,
             "CandidateGenes"=CandidateGenes)
  # WM <-list(WM=Multiplex_Heterogeneous_Matrix,
  #           GeneConnectivity=Connectivity[["gene"]], PhenoConnectivity=Connectivity[["pheno"]],
  #           LG=LG, LP=LP, M=M, N=N, Parameters=Parameters,
  #           gene_pool_nodes_sorted=gene_pool_nodes_sorted,
  #           pheno_pool_nodes_sorted=pheno_pool_nodes_sorted,
  #           CandidateGenes=CandidateGenes, WM_ID = WM_ID[["ID"]], WM_Layers=WM_ID[["LayerNames"]])

  #names(WM) <- c("WM", "GeneConnectivity", "PhenoConnectivity", "LG", "LP", "M", "N", "Parameters",
                 #"gene_pool_nodes_sorted", "pheno_pool_nodes_sorted", "CandidateGenes", "WM_ID", "WM_Layers")
  # saveit(WM=Multiplex_Heterogeneous_Matrix,
  #        GeneConnectivity=Connectivity[["gene"]], PhenoConnectivity=Connectivity[["pheno"]],
  #        LG=LG, LP=LP, M=M, N=N, Parameters=Parameters,
  #        gene_pool_nodes_sorted=gene_pool_nodes_sorted,
  #        pheno_pool_nodes_sorted=pheno_pool_nodes_sorted,
  #        CandidateGenes=CandidateGenes,
  #        WM_ID = WM_ID[["ID"]], WM_Layers=WM_ID[["LayerNames"]],
  #        file = paste0("WM_", WM_ID[["ID"]], ".rda"))

  cat("Time to Start-Finish : ", format(Sys.time()-global_t1), "\n")
  # if there are still some paralllel workers stop and force sequential
  registerDoSEQ()
  #return(WM_ID[["ID"]])
  return(WM)
}

assign.group.to.connectivityDF <- function(ConnectivityDF, no.groups){
  chunk.size <- ceiling(nrow(ConnectivityDF)/no.groups)
  groups <- rep(1:no.groups, each = chunk.size, length.out = nrow(ConnectivityDF))
  ConnectivityDF$Group <- groups
  ConnectivityDF
}
generate.random.seeds <- function(Seeds, ConnectivityDF, S=1000, no.groups=10, replace=FALSE){
  seed.set.size <- length(Seeds)
  sample_size <- ceiling((S/no.groups)*seed.set.size)
  set.seed(1)
  # Stratified Sample "sample_size" nodes from each group as Random Seeds
  RandomSeeds <- ConnectivityDF[which(!ConnectivityDF$Node %in% Seeds),] %>% group_by(Group) %>% sample_n(sample_size, replace = replace)

  # We are creating "seed.set.size" length seed sets by taking "nodes" from each group
  # To this end, we determine a split order and sort the RandomSeeds DF wrt to this "order" column
  order_vec <- vector()
  for(i in 1:no.groups){
    order_vec <- c(order_vec, seq(from=i, to=S*seed.set.size, by=no.groups))
  }
  RandomSeeds$Order <- as.factor(order_vec)
  RandomSeeds <- RandomSeeds[ order(RandomSeeds$Order), ]

  # split the sorted DF into "seed.set.size" length vectors
  RandomSeeds <- split(RandomSeeds$Node, (seq(nrow(RandomSeeds))-1) %/% seed.set.size)

  RandomSeeds
}

#' stratified sampling by using connectivity of the matrix
#'  1. Sort genes by their degree
#'  2. Divided the Genes into 10 groups
#'  3. Get the groups of SeedGenes
#'  4. Random sample from genes from the same degree of genes
#'
#' @param output_dir File directory to output the data.
#' @param WM_ID Id of the walk matrix used.
#' @param ResultID
#' @param GeneSeeds
#' @param GeneConnectivityDF
#' @param PhenoSeeds
#' @param PhenoConnectivityDF
#' @param S
#' @param no.groups.gene
#' @param no.groups.pheno
#'
#' @return
#' @export
#'
#' @examples
generate.random.seed.vector <- function(output_dir, WM_ID, GeneSeeds, GeneConnectivityDF, PhenoSeeds,
                                        PhenoConnectivityDF, S=1000, no.groups.gene=10, no.groups.pheno=5){

  if(length(GeneSeeds)==0 && length(PhenoSeeds)==0) stop("No seeds provided!")
  GeneConnectivity <- assign.group.to.connectivityDF(GeneConnectivity, no.groups=no.groups.gene)
  PhenoConnectivity <- assign.group.to.connectivityDF(PhenoConnectivity, no.groups=no.groups.pheno)


  # sink(paste0(output_dir, "/WM_", WM_ID, "_Connectivity.txt"), append=FALSE)
  # print("Summary of Degree Connectivity of Genes by Group:")
  # print(GeneConnectivity %>%group_by(Group) %>% summarise(n=n(), min=min(Degree), max=max(Degree), mean = mean(Degree)))
  # print("------")
  # print("Summary of Degree Connectivity of Phenoivars by Group:")
  # print(PhenoConnectivity %>%
  #   group_by(Group) %>%
  #   summarise(n=n(), min=min(Degree), max=max(Degree), mean = mean(Degree)))
  # sink()

  #RandomSeeds <- list()
  RandomGeneSeeds <- list()
  if(length(GeneSeeds)!=0){
    RandomGeneSeeds <- generate.random.seeds(GeneSeeds, GeneConnectivity, S, no.groups.gene, TRUE)
    if(length(GeneSeeds)!=1 && any(duplicated(RandomGeneSeeds[1:S]))) warning("WARN: Duplicated random 'gene' seeds generated!")
    #RandomSeeds <- RandomGeneSeeds
  }

  RandomPhenoSeeds <- list()
  if(length(PhenoSeeds)!=0){
    RandomPhenoSeeds <- generate.random.seeds(PhenoSeeds, PhenoConnectivity, S, no.groups.pheno, TRUE)
    if(length(PhenoSeeds)!=1 && any(duplicated(RandomPhenoSeeds[1:S]))) warning("WARN: Duplicated random 'phenotype' seeds generated!")
    #RandomSeeds <- RandomGeneSeeds
  }

  # if(length(GeneSeeds)!=0 & length(PhenoSeeds)!=0){
  #   RandomSeeds <- mapply(c, RandomGeneSeeds, RandomPhenoSeeds, SIMPLIFY=FALSE)
  # }else if(length(RandomSeeds)==0) stop("No seeds provided!")

  return(list(gene=RandomGeneSeeds,pheno=RandomPhenoSeeds))

}

#' Generate S random seeds as length of GeneSeeds excluding  actual Seeds
#'  Run the algorithm for N times, store N result vector
#'   If (the rank of the gene in random seed <= (actual seed rank + Offset) ) then
#'     it will be counted as 1, otherwise 0
#'
#' @param RWGeneRanks
#' @param Rand_Seed_Gene_Rank
#' @param S
#' @param no.cores
#'
#' @return
#' @export
#'
#' @examples
calculatePvalues <- function(RWGeneRanks, Rand_Seed_Gene_Rank, output_dir, no.cores=15){
  #t <- Sys.time()
  S <- ncol(Rand_Seed_Gene_Rank)/3
  cl <- makeCluster(no.cores)
  registerDoParallel(cl)
  dfRanks <- foreach(i = 1:(S)) %dopar% {
    # traverse all gene names and get their ranks in random run result
    rand_ranks <- sapply(RWGeneRanks$Gene, function(gene){ which(Rand_Seed_Gene_Rank[,2+3*(i-1)] %in% gene) })

    # if genes are in the random seed set for this run then there are no rank for them
    idx <- !(sapply(rand_ranks, length))
    rand_ranks[idx] <- NA
    dfRank <- unname(unlist(rand_ranks))
    dfRank
  }
  stopCluster(cl)
  #getDoParWorkers()
  # create dfRanks DF from list output of random seed ranks
  dfRanks <- as.data.frame(bind_cols(dfRanks))

  # add the Gene names as the first column
  dfRanks <- cbind(RWGeneRanks$Gene,RWGeneRanks$Score, dfRanks, stringsAsFactors=FALSE)
  colnames(dfRanks)[1:2] <- c("Gene", "Score")

  rank.offset <- seq(10, 100, by=10)
  # create offset rank columns by adding offset vector for comparison
  dfRanks <- cbind(dfRanks, replicate(length(rank.offset), as.numeric(rownames(dfRanks))) + t(replicate(nrow(dfRanks), rank.offset)))
  colnames(dfRanks)[(S+3):(S+length(rank.offset)+2)] <- paste0("Rank", rank.offset)

  # calculate P values by comparing (random seed rank + offset value) vs (actual seed rank)
  for(i in 1:length(rank.offset)){
    dfRanks$P <- base::rowMeans(dfRanks[, 3:(S+2)] < (dfRanks[,S+2+i]), na.rm = TRUE)

    # rename column P by using offset value e.g. as P10, P20 ...
    dfRanks <- dplyr::rename(dfRanks, (!!as.name(paste0("P",rank.offset[i]))):=P)
  }


  # Calculate Median and Average ranks of genes for Random Seeds
  dfRanks$Med <- apply(dfRanks[,3:(2+S)], 1, median)
  dfRanks$Ave <- rowMeans(dfRanks[,3:(2+S)])
  dfRanks
}

randomWalkRestartsBatch <- function(Walk_Matrix, GeneSeedsList, PhenoSeedsList, N, LG, LP, eta, tau, phi, r, funcs, no.cores=4){
  # t <- Sys.time()
  cl <- makeCluster(no.cores)
  registerDoParallel(cl)
  seedsLength <- ifelse(length(GeneSeedsList)!=0, length(GeneSeedsList), length(PhenoSeedsList))
  #funcs <- c('get.seed.scores', 'randomWalkRestarts', 'rank_genes')

  Rand_Seed_Gene_Rank <- foreach (i=1:seedsLength,.combine=cbind,.export=funcs,.packages=c('Matrix')) %dopar% {
    if (length(GeneSeedsList)!=0 && length(PhenoSeedsList)!=0){
      Seeds_Score <- get.seed.scores(GeneSeedsList[[i]], PhenoSeedsList[[i]], eta, LG, LP, tau, phi )
    }else if(length(GeneSeedsList)!=0){
      Seeds_Score <- get.seed.scores(GeneSeedsList[[i]], vector(), eta, LG, LP, tau, phi )
    }else{
      Seeds_Score <- get.seed.scores(vector(), PhenoSeedsList[[i]], eta, LG, LP, tau, phi )
    }

    Rand_Seed_Res <- randomWalkRestarts(Walk_Matrix, r, Seeds_Score)
    Rand_Seed_Gene_Rank <- rank_genes(N, LG, Rand_Seed_Res, ifelse(length(GeneSeedsList)!=0, GeneSeedsList[[i]], vector()))

    return(Rand_Seed_Gene_Rank)
  }
  stopCluster(cl)
  #cat("Time to Random Walk Batch  : ", format(Sys.time()-t), "\n")

  return(Rand_Seed_Gene_Rank)
}
