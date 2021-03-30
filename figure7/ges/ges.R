library(pcalg)
library(graph)

#do joint estimation given single data
ges.alg <- function(data, k){
	#set up the scoring function
#	score <- new("GaussL0penObsScore", data=data, lambda = 0.5*log(nrow(data)), intercept = FALSE, use.cpp = TRUE)
#	p <- ncol(data)
#
#	#estimate essential graph
#	essgraph <- new("EssGraph", nodes=sapply(1:p, toString), score=score)
#	#start the forward phase
#	essgraph$caus.inf(algorithm="GIES-F", maxSteps=k)
#	#start the backward phase
#	essgraph$caus.inf(algorithm="GIES-B", maxSteps=k)
#	#start the turning phase
#	essgraph$caus.inf(algorithm="GIES-T", maxSteps=k)
#	return(essgraph)

	score <- new("GaussL0penObsScore", data=data, intercept = FALSE, lambda = k * log(ncol(data)), use.cpp = TRUE)

	#estimate essential graph
	ges.fit <- ges(score)
	return(ges.fit$essgraph)
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

#get the estiated graph
#alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7)
#alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
#alphas <- c(0, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 170, 200, 225, 250, 275, 300)
alphas <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 2.2, 2.4, 2.6, 2.8, 3, 4, 5, 6, 7, 8, 16)
grspdag.list <- lapply(alphas, function(alpha) ges.alg(data.list[[1]], alpha))
grspdag.list <- lapply(grspdag.list, function(grspdag) {tt <- as(grspdag, "matrix"); tt[tt == t(tt)] <- 0; as(tt, "graphNEL")})
save(alphas, grspdag.list, file=paste("result/", basename(args[4]), ".rda", sep=""))
