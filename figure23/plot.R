lambdas <- c(0.1, 0.01, 0.001)
neigsize <- c(0.2, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 6, 8) 
#consts <- c(0, 0.25, 0.5, 0.75)
consts <- c(0.25)
cols <- rainbow(3)
ltys <- c(1, 2, 4, 3)
n <- 100000

for(const in consts){
	geometry <- read.table(paste("geometry/Results_10nodes_deg_weight_", substr(toString(const * 100 + 100), 2, 3), sep=""), header=TRUE)
	for(lambda in lambdas){
		prefix <- paste("result/lambda_", toString(lambda), "_weight_", toString(const), sep="")
		load(paste(prefix, "_neig_1_methods.rda", sep=""))
		final.curve <- curve
		for(i in 2:3){
			load(paste(prefix, "_neig_", toString(i), "_methods.rda", sep=""))
			final.curve <- append(final.curve, curve)
		}
		curve <- final.curve
		other.curve <- curve

		load(paste(prefix, "_neig_1_sp.rda", sep=""))
		final.curve <- curve
		for(i in 2:3){
			load(paste(prefix, "_neig_", toString(i), "_sp.rda", sep=""))
			final.curve <- append(final.curve, curve)
		}
		sp.curve <- final.curve

		load(paste("result/weight_", toString(const), "_neig_1_n_", toString(n), "_ges.rda", sep=""))
		final.curve <- curve
		for(i in 2:3){
			load(paste("result/weight_", toString(const), "_neig_", toString(i), "_n_", toString(n), "_ges.rda", sep=""))
			final.curve <- append(final.curve, curve)
		}
		ges.curve <- final.curve

		load(paste(prefix, "_neig_1_sp_d4.rda", sep=""))
		final.curve <- curve
		for(i in 2:3){
			load(paste(prefix, "_neig_", toString(i), "_sp_d4.rda", sep=""))
			final.curve <- append(final.curve, curve)
		}
		sp.d4.curve <- final.curve

		load(paste(prefix, "_neig_1_sp_d5.rda", sep=""))
		final.curve <- curve
		for(i in 2:3){
			load(paste(prefix, "_neig_", toString(i), "_sp_d5.rda", sep=""))
			final.curve <- append(final.curve, curve)
		}
		sp.d5.curve <- final.curve

		load(paste(prefix, "_neig_1_sp_infinity.rda", sep=""))
		final.curve <- curve
		for(i in 2:3){
			load(paste(prefix, "_neig_", toString(i), "_sp_infinity.rda", sep=""))
			final.curve <- append(final.curve, curve)
		}
		sp.infinity.curve <- final.curve

		#proportion of consistent DAGs
		prefix <- paste("figure/lambda_", toString(lambda), "_weight_", toString(const), sep="")
		png(paste(prefix, "_consist.png", sep=""))
		par(mar=c(5,5.5,3,1))
		consis <- sapply(other.curve, function(x) x$prop)
		plot(x=neigsize, y=consis[1,], ylim=c(0, 1), axes=FALSE, xlab="Expected neighbourhood size", ylab="Proportion of consistent simulations", 
			col=cols[1], type="o", cex.lab=2, lwd=3, pch=NA, lty=ltys[1])
		axis(1, las=1, at=0:4 * 2, cex.axis=1.4)
		axis(2, las=1, at=0:5 / 5, cex.axis=1.4)
		box()
		#PC as background
		lines(x=neigsize, y=consis[2,], col="gray", type="o", lwd=3, pch=NA, lty=2)
		#many options for greedy SP
		consis <- sapply(sp.curve, function(x) x$prop)
		d4.consis <- sapply(sp.d4.curve, function(x) x$prop)
		d5.consis <- sapply(sp.d5.curve, function(x) x$prop)
		dinfinity.consis <- sapply(sp.infinity.curve, function(x) x$prop)
		for (i in 1:3){
			lines(x=neigsize, y=d5.consis[i,], col=cols[i], type="o", lwd=3, pch=NA, lty=ltys[2])
			lines(x=neigsize, y=d4.consis[i,], col=cols[i], type="o", lwd=3, pch=NA, lty=ltys[3])
			lines(x=neigsize, y=consis[i*2-1,], col=cols[i], type="o", lwd=3, pch=NA, lty=ltys[4])
		}
		for (i in 2:3){
			lines(x=neigsize, y=dinfinity.consis[i,], col=cols[i], type="o", lwd=3, pch=NA, lty=ltys[1])
		}
		#GES as background
		consis <- sapply(ges.curve, function(x) x$prop)
		lines(x=neigsize, y=consis, col="gray", type="o", lwd=3, pch=NA, lty=3)
		#the hard constraint line
		lines(x=neigsize, y=(1 - geometry$prop[geometry$lambda == lambda]), col="black", type="o", lwd=3, lty=2, pch=NA)
		lname <- c("GES")
		lname <- c(lname, sapply(c(10,5,1), function(i) sapply(c("infinity",5,4,1), function(j) paste("greedy SP, r=", toString(i), ", d=", toString(j), sep=""))))
		lname <- c(lname, "PC")
		#if(lambda == 0.1)
		#	legend("topright", legend=c(lname, "proportion of faithful DAGs"), col=c(cols, "black"), lty=c(rep(1,7),2), lwd=rep(3,8))
		dev.off()
		if(lambda == 0.1){
			png("figure/lengend.png")
			par(mar=c(5,5.5,3,1))
			plot(x=0, y=0)
			legend("topleft", legend=c(lname, "proportion of faithful DAGs"), col=c("gray", rep(cols[3], 4), rep(cols[2], 4), rep(cols[1], 4), "gray", "black"), 
				lty=c(3, ltys, ltys, ltys, 2, 2), lwd=rep(2.4,15), cex=1.5)
			dev.off()
		}

		#average SHD value
		png(paste(prefix, "_shd.png", sep=""))
		par(mar=c(5,5.5,3,1))
		other.aveshd <- sapply(other.curve, function(x) x$aveshd)
		aveshd <- sapply(sp.curve, function(x) x$aveshd)
		d4.aveshd <- sapply(sp.d4.curve, function(x) x$aveshd)
		d5.aveshd <- sapply(sp.d5.curve, function(x) x$aveshd)
		dinfinity.aveshd <- sapply(sp.infinity.curve, function(x) x$aveshd)
		ges.aveshd <- sapply(ges.curve, function(x) x$aveshd)
		plot(x=neigsize, y=other.aveshd[1,], ylim=c(0, 35), axes=FALSE, xlab="Expected neighbourhood size", ylab="Average SHD", 
			col=cols[1], type="o", cex.lab=2, lwd=3, pch=NA, lty=ltys[1])
		axis(1, las=1, at=0:4 * 2, cex.axis=1.4)
		axis(2, las=1, at=0:7 * 5, cex.axis=1.4)
		box()
		lines(x=neigsize, y=other.aveshd[2,], col="gray", type="o", lwd=3, pch=NA, lty=2)
		for (i in 1:3){
			lines(x=neigsize, y=d5.aveshd[i,], col=cols[i], type="o", lwd=3, pch=NA, lty=ltys[2])
			lines(x=neigsize, y=d4.aveshd[i,], col=cols[i], type="o", lwd=3, pch=NA, lty=ltys[3])
			lines(x=neigsize, y=aveshd[2*i-1,], col=cols[i], type="o", lwd=3, pch=NA, lty=ltys[4])
		}
		lines(x=neigsize, y=ges.aveshd, col="gray", type="o", lwd=3, pch=NA, lty=3)
		for (i in 2:3)
			lines(x=neigsize, y=dinfinity.aveshd[i,], col=cols[i], type="o", lwd=3, pch=NA, lty=ltys[1])
		lname <- c("GES")
		lname <- c(lname, sapply(c(1,5,10), function(i) sapply(c(1,4,5,"infinity"), function(j) paste("greedy SP, r=", toString(i), ", d=", toString(j), sep=""))))
		lname <- c(lname, "PC")
		if(lambda == 0.001)
			legend("topleft", legend=lname, col=c("gray", rep(cols[1], 4), rep(cols[2], 4), rep(cols[3], 4), "gray"), lty=c(3, rev(ltys), rev(ltys), rev(ltys), 2), lwd=rep(2,14), cex=1)
		dev.off()
	}
}
