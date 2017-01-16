setwd("~/Dropbox/pauc")
library(pROC)

### data for high pAUC marker
set.seed(234)

ncases <- ncontrols <- 10000 #in each group
# generate a vector of outcomes
y <- c(rep(1, ncases), rep(0, ncontrols))
x1offset <- c(rbeta(ncases, .1, .1), rbeta(ncontrols, .5, 100))
x1 <- scale(x1offset + rnorm(ncases+ncontrols, sd=.3))

rocfit <- roc(response=y, predictor=x1)
auc1 <- rocfit$auc
pauc1 <- auc(rocfit, partial.auc=c(1,.8))
auc1
pauc1

smoothed1 <- smooth(rocfit, method="density")

### data for low pAUC marker
set.seed(125)

ncases <- ncontrols <- 10000 #in each group
# generate a vector of outcomes
y <- c(rep(1, ncases), rep(0, ncontrols))
x2offset <- c(rbeta(ncases, 1.5, .3), rbeta(ncontrols, .4, .5))
x2 <- scale(x2offset + rnorm(ncases+ncontrols, sd=.3))

rocfit <- roc(response=y, predictor=x2)
auc2 <- rocfit$auc
roundval <- round(auc2, 2)

pauc2 <- auc(rocfit, partial.auc=c(1,.8))
auc2
pauc2

smoothed2 <- smooth(rocfit, method="density")

pdf("figures/figure1.pdf", width=14, height=7)
par(mfrow=c(1,2))

plot(0, 0, type="n", xlab="1 - Specificity", ylab="Sensitivity", main="Score A", xlim=c(0,1), ylim=c(0,1))
polygon(list(x=c(0,.2,.2,0), y=c(0,0,1,1)), col="lightgray", lwd=.5)
lines(1-smoothed1$sp, smoothed1$se, col="darkorange2", lwd=2)
t1ind <- which.min(abs(smoothed1$sp-.8))
polyvec <- list(x=c(1-smoothed1$sp[length(smoothed1$sp):t1ind], .2),
				y=c(smoothed1$se[length(smoothed1$se):t1ind], 0))
polygon(polyvec, col="darkblue", density=15, angle=125)
abline(0,1,lty=2,col="darkgray")
text(x=.005,y=-.02,labels=expression(italic(t)[0]))
text(x=.205,y=-.02,labels=expression(italic(t)[1]))
text(x=.9,y=.05,labels=paste("AUC =", round(auc1, 2)))
text(x=.9,y=0,labels=substitute(paste('pAUC'['(0,0.2)'], ' = ', pauc1), list(pauc1=format(pauc1, nsmall=2, digits=2))))

plot(0, 0, type="n", xlab="1 - Specificity", ylab="Sensitivity", main="Score B", xlim=c(0,1), ylim=c(0,1))
polygon(list(x=c(0,.2,.2,0), y=c(0,0,1,1)), col="lightgray", lwd=.5)
lines(1-smoothed2$sp, smoothed2$se, col="darkorange2", lwd=2)
t1ind <- which.min(abs(smoothed2$sp-.8))
polyvec <- list(x=c(1-smoothed2$sp[length(smoothed2$sp):t1ind], .2),
				y=c(smoothed2$se[length(smoothed2$se):t1ind], 0))
polygon(polyvec, col="darkblue", density=15, angle=125)
abline(0,1,lty=2,col="darkgray")
text(x=.005,y=-.02,labels=expression(italic(t)[0]))
text(x=.205,y=-.02,labels=expression(italic(t)[1]))
text(x=.9,y=.05,labels=paste("AUC =", round(auc2, 2)))
text(x=.9,y=0,labels=substitute(paste('pAUC'['(0,0.2)'], ' = ', pauc2), list(pauc2=round(pauc2,2))))

dev.off()

