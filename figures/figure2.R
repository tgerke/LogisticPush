# load and inspect simulation results for logistic push
setwd("~/Dropbox/pauc")

# load results from logistic push
load("code/data/lpSim.RData")
lpselect <- out$selectResults
lppauc <- out$paucResults
lpweights <- out$optWeights

rm(out)

# load results from PAUC-GBS
load("code/data/paucgbsSim.RData")
gbsselect <- out$selectResults
gbspauc <- out$paucResults

rm(out)

# load results from standard lasso
load("code/data/l1Sim.RData")
l1select <- out$selectResults
l1pauc <- out$paucResults

rm(out)

summary(l1pauc)
sd(l1pauc)
t.test(l1pauc)
summary(lppauc)
sd(lppauc)
t.test(lppauc)
summary(gbspauc)
sd(gbspauc)
t.test(gbspauc)

lpselect[1:6]
gbsselect[1:6]
l1select[1:6]
mean(lpselect[7:506])
mean(gbsselect[7:506])
mean(l1select[7:506])

pdf("figures/selectionrates.pdf", height=7, width=11)
par(oma=c(0,0,2,0), xpd=NA)
barplot(rbind(lpselect[1:6], gbsselect[1:6], l1select[1:6])/100, beside=TRUE, ylim=c(0,1), 
		xlab="Marker", ylab="Selection rates", 
		names.arg=c("", "X1", "", "", "X2", "", "", "X3", "", "", "X4", "", "", "X5", "", "", "X6", ""), 
		col=gray.colors(3))
lines(x=c(12.5, 12.5), y=c(-.05,1.05), lty=2)		
text(6.5, .95, "Specific markers")
text(18.5, .95, "Non-specific markers")
legend(par("usr")[2]-4, 1.2, legend=c("Logistic push", "pAUC-GBS", "Lasso"), fill=gray.colors(3))		
dev.off()

summary(lpweights)
plot(hist(lpweights))