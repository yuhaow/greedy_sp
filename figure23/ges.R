library(pcalg)
library(graph)

set.seed(1)

#set up base parameters
simulation <- function(dag.list, data.list){

	#simulation of single random DAG
	dag.alg <- function(data, dag){
		#GES
		score <- new("GaussL0penObsScore", data=data, intercept = FALSE, use.cpp = TRUE)
		gesdag <- as(ges(score)$essgraph, "graphNEL")
		#SHD value
		a <- c(shd(dag$cpdag, gesdag))
		#return the proportion of consistent estimation
		cpdag <- as(dag$cpdag, "matrix")
		b <- c(sum(cpdag != as(gesdag, "matrix")) == 0)
		return(c(b,a))
	}

	#get table of faithful DAGs
	dagnum <- length(dag.list)
	result <- sapply(1:dagnum, function(j) dag.alg(data.list[[j]], dag.list[[j]]))
	faithprop <- sum(result[1,]) / dagnum
	aveshd <- sum(result[2,]) / dagnum
	return(list(prop=faithprop, aveshd=aveshd, result=result))
}

#draw pictures
args <- commandArgs()
load(args[4])

curve <- lapply(1:length(dag.list), function(i) simulation(dag.list[[i]], data.list[[i]]))
suffix <- strsplit(basename(args[4]), ".rda")[[1]][1]
save(curve, file=paste("result/", suffix, "_ges.rda", sep=""))
