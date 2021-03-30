library(graph)
library(pcalg)

#get true positive rate
tprate <- function(est.dag, true.dag){
	essmat <- as(est.dag, "matrix")
	edgenum <- sum(essmat | t(essmat)) / 2
	dag <- as(true.dag, "matrix")
	tp <- sum(essmat & dag)
	skel.tp <- sum((essmat | t(essmat)) & dag)
	return(list(tp=tp, fp=edgenum-tp, skel.tp=skel.tp, skel.fp=edgenum-skel.tp))
}

#plot new figures
get.tpplot <- function(tp.plot, prefix, title){
	maxy = max(sapply(tp.plot, function(t) max(t$tp)))
	maxx = max(sapply(tp.plot, function(t) max(t$fp)))
	png(paste("figure/", prefix, "_", title, "_moralize.png", sep=""))
	par(mar=c(5,5.5,3,1))
	plot(0,0, xlim=c(0, maxx), ylim=c(0, maxy), axes=FALSE,  xlab="Number of false positives", ylab="Number of true positives", cex.lab=2, type="n", cex=2)
	axis(1, las=1, cex.axis=1.4)
	axis(2, las=1, cex.axis=1.4)
	box()
	for(i in 1:length(tp.plot))
		lines(x=tp.plot[[i]]$fp, y=tp.plot[[i]]$tp, col=cols[i], type="p", pch=pchs[i], cex=2)
	legend("bottomright", legend=c("high-dim PC on moral graph", "high-dim GES on moral graph", "high-dim greedy SP on moral graph"), col=cols, pch=pchs, cex=1.5)
	dev.off()
}

ps <- c(100)
ns <- c(300)
#neighs <- c(0.2, 1, 2)
neighs <- c(0.2)
cols <- c("blue", "red", "green3")
pchs <- c(1, 20, 3)

for(p in ps) for(n in ns) for(neigh in neighs){

	#get the prefix of data
	prefix <- paste("p_", toString(p), "_neigh_", toString(neigh), "_n_", toString(n), ".rda", sep="")
	print(prefix)

	#set up the significance levels
	#alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7) # neigh = 1, 2
	alphas <- c(0.0001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7) # neigh = 0.2
	#get results for PC moralize
	pc.moral.list <- lapply(alphas, function(alpha){
		load(paste("pc_moralize/result/", prefix, "_alpha_", toString(alpha), ".rda", sep=""))
		a <- sapply(1:length(dag.list), function(i) tprate(pcdag.list[[i]], dag.list[[i]]))
		list(tp=mean(unlist(a["tp",])), fp=mean(unlist(a["fp",])), skel.tp=mean(unlist(a["skel.tp",])), skel.fp=mean(unlist(a["skel.fp",])))
	})
	alphas <- c(0.000001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
	#get results for SP moralize
	sp.moral.list <- lapply(alphas, function(alpha){
		load(paste("sp_moralize/result/", prefix, "_alpha_", toString(alpha), "_restart.rda", sep=""))
		a <- sapply(1:length(dag.list), function(i) tprate(grspdag.list[[i]], dag.list[[i]]))
		list(tp=mean(unlist(a["tp",])), fp=mean(unlist(a["fp",])), skel.tp=mean(unlist(a["skel.tp",])), skel.fp=mean(unlist(a["skel.fp",])))
	})

	#set up the penalization constants
	consts <- c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16)
	#get results for GES moralize
	ges.moral.list <- lapply(consts, function(const){
		load(paste("ges_moralize/result/", prefix, "_const_", toString(const), ".rda", sep=""))
		a <- sapply(1:length(dag.list), function(i) tprate(gesdag.list[[i]], dag.list[[i]]))
		list(tp=mean(unlist(a["tp",])), fp=mean(unlist(a["fp",])), skel.tp=mean(unlist(a["skel.tp",])), skel.fp=mean(unlist(a["skel.fp",])))
	})

	#plot the results, first comes the rates for directed edges
	#di.tpplot <- lapply(list(pc.moral.list, ges.moral.list, sp.moral.list), function(tp.list) 
	#	list(tp=sapply(tp.list, function(t) t$tp), fp=sapply(tp.list, function(t) t$fp)))
	#get.tpplot(di.tpplot, prefix, "directed")

	#plot the results, second comes the rates for skeleton
	skel.tpplot <- lapply(list(pc.moral.list, ges.moral.list, sp.moral.list), function(tp.list) 
		list(tp=sapply(tp.list, function(t) t$skel.tp), fp=sapply(tp.list, function(t) t$skel.fp)))
	get.tpplot(skel.tpplot, prefix, "skeleton")
}
