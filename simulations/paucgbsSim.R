Sys.info()
sessionInfo()

# this simulation will evaluate variable selection rates for pauc-gbs and determine pauc peformance in test data
setwd("~/Dropbox/pauc")
library(mboost)
library(pROC)
source("code/helpers/paucgbsRisk.R")
source("code/helpers/genData.R")

nsim <- 100
upperFPR <- .2
# a matrix which will count the number of times each variable is selected over nsim replications
selectResults <- rep(0, 506)
paucResults <- rep(0, nsim)

for (i in 1:nsim) {
	simData <- 	genData(50, 50, 3, 3, 500, seed=i*2)
	simData$y <- as.factor(simData$y)
	sigma <- min(table(simData$y))^(-1/4)
	INT <- rep(1, nrow(simData))

	# specify the base-learners for component-wise gradient boosting
	formula1 <- as.formula(paste("y ~", paste("bols(", names(simData)[-1], ",
			intercept = FALSE)", sep = "",  collapse = " + ")))
			
	# run component-wise gradient boosting with FPR range [0,0.2]
	model1 <- gamboost(formula1, data = simData, family = PAUC(fprup = upperFPR,
			sigma = sigma), control=boost_control(trace = FALSE,
			mstop = 50, nu = 0.1))
	# generate subsamples for stability selection
	cv5f1 <- cv(model.weights(model1), type = "kfold", B = 5)
	cv5f2 <- cv(model.weights(model1), type = "kfold", B = 5)
	cv5f3 <- cv(model.weights(model1), type = "kfold", B = 5)
	cv5f4 <- cv(model.weights(model1), type = "kfold", B = 5)
	cv5f <- cbind(cv5f1, cv5f2, cv5f3, cv5f4)
	
	# run stability selection
	STAB <- stabsel(model1, FWER = 0.1, cutoff = 0.9, folds = cv5f)
	
	# re-run component-wise gradient boosting, this time using the selected variables only
	if(length(names(STAB$selected)) > 0) {
		blsnames2 <- paste(names(STAB$selected), sep = "", collapse = "+")
		formula2 <- as.formula(paste("y ~ ", paste(blsnames2, sep = "+")))
	} else {
		blsnames2 <- ""
		formula2 <- as.formula(paste("y ~ bols(INT, intercept = FALSE)"))
	}
	model1 <- gamboost(formula2, data=simData, family=PAUC(fprup = upperFPR,
		sigma = sigma), control = boost_control(trace = FALSE, mstop = 200, nu = 0.1))
	selectinds <- which(sapply(names(simData)[-1], function(x) grep(x, names(STAB$selected)))==TRUE)
	selectResults[selectinds] <- selectResults[selectinds]+1

	simData2 <- genData(50, 50, 3, 3, 500, seed=i*200)
	preds <- predict(model1, newdata=simData2)
	paucResults[i] <- auc(roc(response=simData2$y, predictor=preds), partial.auc=c(1,1-upperFPR))
}

out <- list(selectResults=selectResults, paucResults=paucResults)
save(out, file="data/paucgbsSim.RData")