library(pcalg)
library(graph)
library(MASS)

#parameter used for generating data
ps <- c(8)
ns <- c(1000, 10000)
neighs <- c(4)
dagnum <- 100

#set random seed
set.seed(1)

#generate random graphs
for(p in ps){
for(neigh in neighs){
	#function for generating random dags
	rnddag <- function(){
		dag <- randomDAG(p, prob = neigh / (p - 1), lB = 0.25, uB = 1)
		dag <- as(dag, "matrix") * matrix(sample(c(-1, 1), p * p, replace=TRUE), nrow=p)
		dag <- lapply(1:p, function(i) list(edges=which(dag[i,] != 0), weights=dag[i, dag[i,] != 0]))
		names(dag) <- 1:p
		dag <- new("graphNEL", node=sapply(1:p, toString), edgeL=dag, edgemode="directed")
		return(dag)
	}

	#get list of dags
	dag.list <- lapply(1:dagnum, function(i) rnddag())
	for(n in ns){
		data.list <- lapply(dag.list, function(dag) mvrnorm(n, mu=rep(0,p), Sigma=trueCov(dag)))
		title <- paste("data/p_", toString(p), "_neigh_", toString(neigh), "_n_", toString(n), ".rda", sep="")
		save(data.list, dag.list, dagnum, title, file=title)
	}
}}
