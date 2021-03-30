library(lattice)

source("lines.R")

p = 8
n = 10000	
min_weight = 0.25
alpha = c(0.01, 0.001, 0.0001)

ntests = 6
sparsity = c(0.2, 1:(p-1))

# read in file Plot_Figs_points.rtf

data=rep(0,length(sparsity)*length(alpha)*ntests)
dim(data)=c(length(sparsity),length(alpha),ntests)
dimnames(data)=list(sparsity, alpha, c("SP", "greedy SP", "PC", "SGS", "MMHC", "GES"))

data[,1,1] = Results_SP_8_10000_0.01
data[,2,1] = Results_SP_8_10000_0.001
data[,3,1] = Results_SP_8_10000_0.0001
data[,1,2] = Results_PC_8_10000_0.01
data[,2,2] = Results_PC_8_10000_0.001
data[,3,2] = Results_PC_8_10000_0.0001
data[,1,3] = Results_SGS_8_10000_0.01
data[,2,3] = Results_SGS_8_10000_0.001
data[,3,3] = Results_SGS_8_10000_0.0001
data[,1,6] = Results_GSP_8_10000_0.01
data[,2,6] = Results_GSP_8_10000_0.001
data[,3,6] = Results_GSP_8_10000_0.0001
for (i in 1:3){
	 data[,i,4] = Results_MMHC_8_10000
	 data[,i,5] = Results_GES_8_10000
}

### To get line plots for paper

for(j in 1:3){

	png(paste("figures/n10000_alpha_", alpha[j], ".png", sep=""))

	par(mar=c(5,5.5,3,1))
	plot(sparsity, data[,j,1], xlim=c(-0.05,(p-0.95)), ylim=c(0,1), axes=FALSE, xlab="Expected neighborhood size", ylab="Proportion of simulations", cex.lab=2, line=3.5, pch=NA)

	axis(1, las=1, at=0:p, cex.axis=1.4)
	axis(2, las=1, at=c(0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), cex.axis=1.4)

	lines(c(0,0), c(-0.1,1.1), lwd=0.4, col="grey")
	lines(c(1,1), c(-0.1,1.1), lwd=0.4, col="grey")
	lines(c(2,2), c(-0.1,1.1), lwd=0.4, col="grey")
	lines(c(3,3), c(-0.1,1.1), lwd=0.4, col="grey")
	lines(c(4,4), c(-0.1,1.1), lwd=0.4, col="grey")
	lines(c(5,5), c(-0.1,0.673), lwd=0.4, col="grey")
	lines(c(5,5), c(1.01,1.1), lwd=0.4, col="grey")
	lines(c(6,6), c(-0.1,0.673), lwd=0.4, col="grey")
	lines(c(6,6), c(1.01,1.1), lwd=0.4, col="grey")
	lines(c(7,7), c(-0.1,0.673), lwd=0.4, col="grey")
	lines(c(7,7), c(1.01,1.1), lwd=0.4, col="grey")

	box()

	lines(sparsity, data[,j,1], type="o", pch=21, lty=1, col="red", lwd=3)
	lines(sparsity, data[,j,2], type="o", pch=22, lty=3, col="blue", lwd=3)
	lines(sparsity, data[,j,3], type="o", pch=23, lty=5, col="green3", lwd=3)
	lines(sparsity, data[,j,4], type="o", pch=24, lty=4, col="grey25", lwd=3)
	lines(sparsity, data[,j,5], type="o", pch=25, lty=6, col="maroon4", lwd=3)
	lines(sparsity, data[,j,6], type="o", pch=8, lty=2, col="orange", lwd=3)

	legend("topright", c("SP", "greedy SP", "PC","SGS", "GES", "MMHC"), cex=1.4, 
		col=c("red", "orange", "blue", "green3", "maroon4", "grey25"), pch=c(21,8,22,23,25,24), lty=c(1,2,3,5,6,4), lwd=c(3,3,3,3,3,3));

	title(main=paste("8 nodes, nsamples 10000, alpha", alpha[j]), font.main=1, cex.main=2)

	dev.off()
}

data[,1,1] = Results_SP_8_1000_0.01
data[,2,1] = Results_SP_8_1000_0.001
data[,3,1] = Results_SP_8_1000_0.0001
data[,1,2] = Results_PC_8_1000_0.01
data[,2,2] = Results_PC_8_1000_0.001
data[,3,2] = Results_PC_8_1000_0.0001
data[,1,3] = Results_SGS_8_1000_0.01
data[,2,3] = Results_SGS_8_1000_0.001
data[,3,3] = Results_SGS_8_1000_0.0001
for (i in 1:3){
	 data[,i,4] = Results_MMHC_8_1000
	 data[,i,5] = Results_GES_8_1000
}
data[,1,6] = Results_GSP_8_1000_0.01
data[,2,6] = Results_GSP_8_1000_0.001
data[,3,6] = Results_GSP_8_1000_0.0001

for(j in 1:3){

	png(paste("figures/n1000_alpha_", alpha[j], ".png", sep=""))

	par(mar=c(5,5.5,3,1))
	plot(sparsity, data[,j,1], xlim=c(-0.05,(p-0.95)), ylim=c(0,1), axes=FALSE, xlab="Expected neighborhood size", ylab="Proportion of simulations", cex.lab=2, line=3.5, pch=NA)

	axis(1, las=1, at=0:p, cex.axis=1.4)
	axis(2, las=1, at=c(0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), cex.axis=1.4)

	lines(c(0,0), c(-0.1,1.1), lwd=0.4, col="grey")
	lines(c(1,1), c(-0.1,1.1), lwd=0.4, col="grey")
	lines(c(2,2), c(-0.1,1.1), lwd=0.4, col="grey")
	lines(c(3,3), c(-0.1,1.1), lwd=0.4, col="grey")
	lines(c(4,4), c(-0.1,1.1), lwd=0.4, col="grey")
	lines(c(5,5), c(-0.1,0.673), lwd=0.4, col="grey")
	lines(c(5,5), c(1.01,1.1), lwd=0.4, col="grey")
	lines(c(6,6), c(-0.1,0.673), lwd=0.4, col="grey")
	lines(c(6,6), c(1.01,1.1), lwd=0.4, col="grey")
	lines(c(7,7), c(-0.1,0.673), lwd=0.4, col="grey")
	lines(c(7,7), c(1.01,1.1), lwd=0.4, col="grey")

	box()

	lines(sparsity, data[,j,1], type="o", pch=21, lty=1, col="red", lwd=3)
	lines(sparsity, data[,j,2], type="o", pch=22, lty=3, col="blue", lwd=3)
	lines(sparsity, data[,j,3], type="o", pch=23, lty=5, col="green3", lwd=3)
	lines(sparsity, data[,j,4], type="o", pch=24, lty=4, col="grey25", lwd=3)
	lines(sparsity, data[,j,5], type="o", pch=25, lty=6, col="maroon4", lwd=3)
	lines(sparsity, data[,j,6], type="o", pch=8, lty=2, col="orange", lwd=3)

	legend("topright", c("SP", "greedy SP", "PC","SGS", "GES", "MMHC"), cex=1.4, 
		col=c("red", "orange", "blue", "green3", "maroon4", "grey25"), pch=c(21,8,22,23,25,24), lty=c(1,2,3,5,6,4), lwd=c(3,3,3,3,3,3));

	title(main=paste("8 nodes, nsamples 1000, alpha", alpha[j]), font.main=1, cex.main=2)

	dev.off()
}
