library(pcalg)
library(graph)
library(MASS)

#source("generate_random_DAG.R")

#parameter used for generating data
ps <- c(8)
ns <- c(10000, 1000)
neighs <- c(0.2, 1, 2, 3, 4, 5, 6, 7)
dagnum <- 100

#set random seed
set.seed(1)

rndDAG <- function (n, prob, lB, uB)
{
    stopifnot(n >= 2, is.numeric(prob), length(prob) == 1, 0 <= prob, prob <= 1,
              is.numeric(lB), is.numeric(uB), lB <= uB)
    V <- as.character(1:n)
    edL <- as.list(V)
    names(edL) <- V
    nmbEdges <- 0
    genWgt <- function(n, min, max) {
      w <- runif(n, min = min, max = max)
      sgn <- sample(c(-1,1), n, replace = TRUE)
      w*sgn
    }
    for (i in seq(length = n - 2)) {
        listSize <- rbinom(1, n - i, prob)
        nmbEdges <- nmbEdges + listSize
        edgeList <- sample(seq(i + 1, n), size = listSize)
        ## weightList <- runif(length(edgeList), min = lB, max = uB)
        weightList <- genWgt(length(edgeList), min = lB, max = uB)
        edL[[i]] <- list(edges = edgeList, weights = weightList)
    }
    if (rbinom(1, 1, prob) == 1) {
        edL[[n - 1]] <- list(edges = n, weights = genWgt(1, min = lB, max = uB))
    }
    else {
        edL[[n - 1]] <- list(edges = integer(0), weights = numeric(0))
    }
    edL[[n]] <- list(edges = integer(0), weights = numeric(0))
    if (nmbEdges > 0) {
        res <- new("graphNEL", nodes = V, edgeL = edL, edgemode = "directed")
    } else {
        res <- new("graphNEL", nodes = V, edgemode = "directed")
    }
    res
}

#function for generating random dags
randgraph <- function(p, neigh){
	dag <- randomDAG(p, prob = neigh / (p - 1), lB = 0.25, uB = 1)
	dag <- as(dag, "matrix") * matrix(sample(c(-1, 1), p * p, replace=TRUE), nrow=p)
	dag <- lapply(1:p, function(i) list(edges=which(dag[i,] != 0), weights=dag[i, dag[i,] != 0]))
	names(dag) <- 1:p
	dag <- new("graphNEL", node=sapply(1:p, toString), edgeL=dag, edgemode="directed")
	return(dag)
}

#generate random graphs
for(p in ps){
for(neigh in neighs){

	#get list of dags
	dag.list <- lapply(1:dagnum, function(i) randgraph(p, neigh))
	for(n in ns){
		data.list <- lapply(dag.list, function(dag) mvrnorm(n, mu=rep(0,p), Sigma=trueCov(dag)))
		title <- paste("data/p_", toString(p), "_neigh_", toString(neigh), "_n_", toString(n), ".rda", sep="")
		save(data.list, dag.list, dagnum, title, file=title)
	}
}}
