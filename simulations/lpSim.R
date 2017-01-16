# this simulation will evaluate variable selection rates for logistic push and determine pauc peformance in test data
setwd("~/Dropbox/pauc")
library(glmnet)
library(pROC)
source("code/helpers/genData.R")
source("code/helpers/logpush.R")

Sys.info()
sessionInfo()

nsim <- 100
upperFPR <- .2
# a matrix which will count the number of times each variable is selected over nsim replications
selectResults <- rep(0, 506)
paucResults <- rep(0, nsim)
# vector of weights over which to search
wvec <- seq(1,20,by=.25)
optWeights <- rep(0, nsim)

for (k in 1:nsim) {
	simData <- 	genData(50, 50, 3, 3, 500, seed=k*2)

	# set up cross-validation sets for the weights (unstratified)
	#foldids <- sample(rep(1:5, length=nrow(simData)))
	#results <- rep(0, length(wvec))

	# set up stratified cross-validation sets for the weights
	y1foldids <- sample(rep(1:5, length=sum(simData$y==1)))
	y0foldids <- sample(rep(1:5, length=sum(simData$y==0)))	
	results <- rep(0, length(wvec))
	
	lpseed <- 111
	for (i in 1:length(wvec)) {
  		for (j in 1:5) {
  			# unstratified cv case
    			#traininds <- which(foldids!=j)
    			#testinds <- which(foldids==j)
    		
    			# stratified case
    			traininds <- c(which(simData$y==1)[which(y1foldids!=j)], which(simData$y==0)[which(y0foldids!=j)])
			testinds <- c(which(simData$y==1)[which(y1foldids==j)], which(simData$y==0)[which(y0foldids==j)])
    
   			lpfit <- logpush(x=as.matrix(simData[traininds,-1]), 
                     y=simData$y[traininds],
                     wts=ifelse(simData$y[traininds]==1, 1, wvec[i]), 
                     nlambda=100, nfolds=5, upperFPR=upperFPR, seed=lpseed)
		    preds <- as.numeric(predict(lpfit$fit, newx=as.matrix(simData[testinds,-1]),
                     s=lpfit$lambdaopt))
    		results[i] <- results[i] +
    				 pROC:::auc(roc(response=simData$y[testinds], predictor=preds), 
          			 partial.auc=c(1,1-upperFPR))
          	lpseed <- lpseed + 111
  		}
	}
	optWeights[k] <- wvec[which.max(results)]

	weightvec <- ifelse(simData$y==1, 1, optWeights[k])
	lpfit <- logpush(x=as.matrix(simData[,-1]), y=simData$y,wts=weightvec, 
                 nlambda=100, nfolds=5, upperFPR=upperFPR)
	nz <- unlist(predict(lpfit$fit, type="nonzero", s=lpfit$lambdaopt))
	selectResults[nz] <- selectResults[nz] + 1

	simData2 <- genData(50, 50, 3, 3, 500, seed=k*200)
	preds <- as.numeric(predict(lpfit$fit, newx=as.matrix(simData2[,-1]), s=lpfit$lambdaopt, type="response"))
	paucResults[k] <- pROC:::auc(roc(response=simData2$y, predictor=preds), partial.auc=c(1,1-upperFPR))
	print(k)
}


save(list(selectResults=selectResults, paucResults=paucResults, optWeights=optWeights), file="data/lpSim.RData")
