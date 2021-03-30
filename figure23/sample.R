library(pcalg)
library(graph)
library(MASS)

set.seed(1)

#set up base parameters
dag.proc <- function(p, neigh, const){

	#randomly generate dag
	dag <- randomDAG(p, prob = neigh / (p - 1), lB = const, uB = 1)
	dagmat <- as(dag, "matrix") * matrix(sample(c(-1, 1), p * p, replace=TRUE), nrow=p)
	dagmat <- (diag(rep(1, p)) - dagmat) %*% t(diag(rep(1, p)) - dagmat)
	#original dag to cpdag
	cpdag <- dag2cpdag(dag)
	return(list(dagmat=dagmat, cpdag=cpdag))
}

#draw pictures
neigsize <- list(c(0.2, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4), c(6), c(8))
consts <- c(0.25)
ns <- c(1000000)
for(const in consts){
	for(neig.idx in 1:length(neigsize)){
		load(paste("data/weight_", toString(const), "_neig_", toString(neig.idx), ".rda", sep=""))
		for(n in ns){
			data.list <- lapply(dag.list, function(dags) lapply(dags, function(dag) mvrnorm(n, mu=rep(0,10), Sigma=solve(dag$dagmat))))
			title <- paste("data/weight_", toString(const), "_neig_", toString(neig.idx), "_n_", toString(n), ".rda", sep="")
			save(data.list, dag.list, title, file=title)
		}
	}
}
