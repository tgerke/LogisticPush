# this function will generate data follows distributions set forth by Schmid et al.
# nSpMarkers is the number of specific markers 
# nMarkers is the number of informative (but not highly specific) markers
# nNoise is the number of noise variables
genData <- function(ncases, ncontrols, nSpMarkers, nMarkers, nNoise, seed) {
	set.seed(seed)

	# generate a vector of outcomes
	y <- c(rep(1, ncases), rep(0, ncontrols))
	
	# generate a matrix of specific marker values
	spMat <- matrix(0, nrow=ncases+ncontrols, ncol=nSpMarkers)
	spOffset <- c(rbeta(ncases, .1, .1), rbeta(ncontrols, .5, 100))
	for (i in 1:ncol(spMat)) {
		spMat[,i] <- scale(spOffset + rnorm(ncases+ncontrols, sd=.3))
	}
	
	# generate matrix of informative marker values
	markMat <- matrix(0, nrow=ncases+ncontrols, ncol=nMarkers)
	markOffset <- c(rbeta(ncases, 1.5, .3), rbeta(ncontrols, .4, .5))
	for (i in 1:ncol(markMat)) {
		markMat[,i] <- scale(markOffset + rnorm(ncases+ncontrols, sd=.3))
	}	
	
	# generate matrix of noise variables
	noiseMat <- matrix(runif(nNoise*(ncases+ncontrols)), ncol=nNoise)
	noiseMat <- apply(noiseMat, 2, scale)
	
	dat <- data.frame(cbind(y, spMat, markMat, noiseMat))
	names(dat) <- c("y", paste("sp", 1:nSpMarkers, sep=""), paste("m", 1:nMarkers, sep=""),
					paste("noise", 1:nNoise, sep=""))
					
	return(dat)
}

# ### for testing
# library(pROC)
# myDat <- genData(50, 50, 3, 3, 500, 8675309)
# rocfit <- roc(response=myDat$y, predictor=myDat$sp1)
# rocfit$auc
# auc(rocfit, partial.auc=c(1,.8))
# par(mfrow=c(1,2))
# plot(rocfit, print.auc=TRUE, 
		 # auc.polygon=TRUE, partial.auc=c(1, 0.8),
         # partial.auc.focus="sp", grid=c(0.1, 0.2), grid.col=c("green", "red"),
         # max.auc.polygon=TRUE, auc.polygon.col="blue", print.thres=FALSE,
         # reuse.auc=FALSE, legacy.axes=FALSE, main="Specific marker", 
         # xlab="True negative rate", ylab="True positive rate")
# rocfit <- roc(response=myDat$y, predictor=myDat$m1)
# rocfit$auc
# auc(rocfit, partial.auc=c(1,.8))
# plot(rocfit, print.auc=TRUE, 
		 # auc.polygon=TRUE, partial.auc=c(1, 0.8),
         # partial.auc.focus="sp", grid=c(0.1, 0.2), grid.col=c("green", "red"),
         # max.auc.polygon=TRUE, auc.polygon.col="blue", print.thres=FALSE,
         # reuse.auc=FALSE, legacy.axes=FALSE, main="Non-specific marker", 
         # xlab="True negative rate", ylab="True positive rate")