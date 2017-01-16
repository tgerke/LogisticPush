# logpush fits a cross-validated logistic push model with the supplied weights
# x is n-by-p matrix of variables
# y is vector of outcome values
logpush <- function(x, y, wts, nlambda, nfolds, upperFPR, seed=111) {
	set.seed(seed)

	# set up stratified foldids for cv
	foldids <- rep(0, length=length(y))
	y1foldids <- sample(rep(1:nfolds, length=sum(y==1)))
	y0foldids <- sample(rep(1:nfolds, length=sum(y==0)))	
	foldids[which(y==1)] <- y1foldids
	foldids[which(y==0)] <- y0foldids
	
	cvfit <- cv.glmnet(x=x, y=y, family="binomial", alpha=1, weights=wts, 
				   nlambda=nlambda, foldid=foldids, keep=TRUE)
	
	cvpauc <- vector("list", length=nfolds)
	for (i in 1:nfolds) {
		inds <- which(cvfit$foldid==i)
		cvpauc[[i]] <- apply(cvfit$fit.preval[,colSums(is.na(cvfit$fit.preval))==0], 2, 
				function(z) {pROC:::auc(roc(response=y[inds], predictor=z[inds]), partial.auc=c(1,1-upperFPR))})
	}
	cvavg <- Reduce("+", cvpauc)/nfolds
	# return FIRST index of max value, which gives the fewest selected variables at that lambda 
	lambdaopt <- cvfit$lambda[which.max(cvavg)]

	return(list(fit=cvfit$glmnet.fit, lambdaopt=lambdaopt, maxpauc=max(cvavg)))	
}

# for testing
#simData <- genData(50, 50, 3, 3, 500, 1234)
#lpfit <- logpush(x=as.matrix(simData[,-1]), y=simData$y, wts=ifelse(simData$y==1, 1, 10), 
#					nlambda=250, nfolds=5, upperFPR=.2)
#predict(lpfit$fit, type="nonzero", s=lpfit$lambdaopt)