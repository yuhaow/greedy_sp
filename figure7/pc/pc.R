library(pcalg)
library(graph)

#input to Greedy SP with random restarts: mat -- precision matrix, order -- initial permutation
#input to Greedy SP with random restarts: mat -- precision matrix, order -- initial permutation
pc.alg <- function(suffstat, alpha){
	p <- nrow(suffstat$C)
	pc.fit <- pc(suffStat=suffstat, indepTest=gaussCItest, alpha=alpha, labels=sapply(1:p, toString), verbose=TRUE)
	return(pc.fit@graph)
}

#get the test significance for estimated edges
conf.mat <- function(dag, order, suffstat){
	p <- ncol(suffstat$C)
	revorder <- sapply(1:p, function(t) which(order==t))
	t <- sapply(1:p, function(j) sapply(1:p, function(i) if(dag[i,j] != 0) 1 - gaussCItest(i, j, order[c(1:(revorder[j]-1))[-revorder[i]]], suffstat) else 0))
	t <- sapply(1:p, function(j) sapply(1:p, function(i) if(dag[i,j] != 0) t[i, j] else t[j, i]))
}

#get data for simulation
args <- commandArgs()
load(args[4])
print(title)

set.seed(1)

#get result for each dag
suffstat <- list(C=cor(data.list[[1]]), n=nrow(data.list[[1]]))

#get the estiated graph
#alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7)
#alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
#alphas <- c(0.75, 0.8, 0.85)
alphas <- c(0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
grspdag.list <- lapply(alphas, function(alpha) pc.alg(suffstat, alpha))
grspdag.list <- lapply(grspdag.list, function(grspdag) {tt <- as(grspdag, "matrix"); tt[tt == t(tt)] <- 0; as(tt, "graphNEL")})
save(alphas, grspdag.list, file=paste("result/", basename(args[4]), ".rda", sep=""))
