#######################################################################################
# this function from PMID 23033406 specifies pAUC as the risk function for boosting

# fprup corresponds to the upper limit of the FPR range
PAUC <- function (fprup = 1, sigma = 0.1) {
    approxGrad <- function(x) {exp(-x/sigma) / (sigma * (1 + exp(-x/sigma))^2)}
    approxLoss <- function(x) {1 / (1 + exp(-x / sigma))}
    
    Family(
        # implement the gradient of PAUC (formulas (2.6) and (2.7))
        ngradient = function(y, f, w = 1) {
        if (!all(w %in% c(0,1)))
            stop(sQuote("weights"), " must be either 0 or 1 for family ", sQuote("PAUC"))
        if (length(w) == 1) w <- rep(1, length(y))
        ind1 <- which(y == 1)
        ind0 <- which(y == -1)
        n1 <- length(ind1)
        n0 <- length(ind0)
        if (length(f) == 1) {f <- rep(f, length(y))}
        f <- f - f[w == 1][1]
        # build weight matrix
        tmp <- matrix(w[ind1], nrow = n0, ncol = n1, byrow = TRUE)
        weightmat <- matrix(w[ind0], nrow = n0, ncol = n1) * tmp
        # differences between "diseased" and "non-diseased"
        M0 <- matrix(-f[ind1], nrow = n0, ncol = n1, byrow = TRUE) + f[ind0]
        M1 <- approxGrad(M0) * weightmat
        M2 <- approxLoss(M0) * weightmat
        denom <- 1 + exp( (colSums(M2) / sum(w[ind0]) - fprup) / sigma )
        ng <- vector(length(y), mode = "numeric")
        ng[ind1] <- colSums(M1) / denom / sigma / (sum(w[ind1]))
        ng[ind0] <- rowSums(- sweep(M1, 2, denom, FUN = "/")) / sigma / sum(w[ind1])
        return(ng)
        },
        
        # implement the smoothed negative PAUC risk (formula (2.5))
        risk = function(y, f, w = 1) {
        if (length(w) == 1) w <- rep(1, length(y))
        ind1 <- which(y == 1)
        ind0 <- which(y == -1)
        n1 <- length(ind1)
        n0 <- length(ind0)
        if (length(f) == 1) {f <- rep(f, length(y))}
        f <- f - f[w == 1][1]
        tmp <- matrix(w[ind1], nrow = n0, ncol = n1, byrow = TRUE)
        weightmat <- matrix(w[ind0], nrow = n0, ncol = n1) * tmp
        M0 <- matrix(-f[ind1], nrow = n0, ncol = n1, byrow = TRUE) + f[ind0]
        M1 <- approxGrad(M0) * weightmat
        M2 <- approxLoss(M0) * weightmat
        num <- 1 + exp(fprup / sigma)
        denom <- 1 + exp( (fprup - colSums(M2) / sum(w[ind0])) / sigma )
        return(-(sum(fprup - sigma * log(num/denom))) / (sum(w[ind1])))
        },
        
        weights = "case",
        
        offset = function(y, w) {0},
        
        check_y = function(y) {
        if (!is.factor(y))
            stop("response is not a factor but ", sQuote("family = PAUCSigma()"))
        if (nlevels(y) != 2)
            stop("response is not a factor at two levels but ", sQuote("family = AUC()"))
        if (length(unique(y)) != 2)
            stop("only one class is present in response.")
        ind1 <- which(y == levels(y)[2])
        ind0 <- which(y == levels(y)[1])
        n1 <- length(ind1)
        n0 <- length(ind0)
        c(-1, 1)[as.integer(y)]
        },
        
        rclass = function(f) (f > 0) + 1,

        name = paste("(1 - Partial AUC)-Loss")
    )
}