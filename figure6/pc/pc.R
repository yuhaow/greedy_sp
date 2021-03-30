library(pcalg)
library(graph)
library(parallel)

set.seed(1)

#PC algorithm 
pc.alg <- function(suffstat){
	p <- nrow(suffstat$C)
	pc.fit <- pc(suffStat=suffstat, indepTest=gaussCItest, alpha=alpha, labels=sapply(1:p, toString), verbose=TRUE)
	return(pc.fit@graph)
}

#get data for simulation
args <- commandArgs()
load(args[4])
print(title)
alpha <- as.numeric(args[5])
#dagnum <- 5

#set up parallels
cor.num <- min(detectCores(), dagnum)
cl <- makeCluster(getOption("cl.cores", cor.num))
clusterExport(cl, c("alpha"))
clusterEvalQ(cl, c(library(pcalg), library(graph)))

#start parallel computation
suffstat.list <- lapply(1:dagnum, function(i) list(C=cor(data.list[[i]]), n=nrow(data.list[[i]])))
pcdag.list <- parLapply(cl, suffstat.list, pc.alg)
stopCluster(cl)

save(pcdag.list , dag.list, file=paste("result/", basename(args[4]), "_alpha_", toString(alpha), ".rda", sep=""))
