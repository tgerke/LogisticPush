setwd("~/Dropbox/pauc/code/LogisticPush")
library(glmnet)
library(pROC)
source("helpers/genData.R")
source("helpers/logpush.R")

set.seed(8675309)

# maximum allowable false positive rate
upperFPR <- .2
# vector of weights over which to search
wvec <- seq(1,20,by=2)

simData <- 	genData(ncases=50, ncontrols=50, nSpMarkers=3, nMarkers=3, 
                    nNoise=500, seed=8675309)
print(simData[c(1:5,96:100), 1:10])

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
    	traininds <- c(which(simData$y==1)[which(y1foldids!=j)], 
                     which(simData$y==0)[which(y0foldids!=j)])
		testinds <- c(which(simData$y==1)[which(y1foldids==j)], 
                     which(simData$y==0)[which(y0foldids==j)])
    
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
   cat("iteration", i, "of", length(wvec), "\n")
}
# CV-averaged results
print(results/5)
# optimal weight
optweight <- wvec[which.max(results)]
print(optweight)

weightvec <- ifelse(simData$y==1, 1, optweight)
lpfit <- logpush(x=as.matrix(simData[,-1]), y=simData$y, wts=weightvec, 
                 nlambda=100, nfolds=5, upperFPR=upperFPR)
coefs <- coef(lpfit$fit, s=lpfit$lambdaopt)
# print selected markers
coefs <- coefs[which(coefs!=0),]
print(coefs)

# simulate new (test) data
simData2 <- genData(50, 50, 3, 3, 500, seed=999)
preds <- as.numeric(predict(lpfit$fit, newx=as.matrix(simData2[,-1]), s=lpfit$lambdaopt, type="response"))
# check validation set performance
pROC:::auc(roc(response=simData2$y, predictor=preds), partial.auc=c(1,1-upperFPR))

