#' Deviance function using the link function value for different families
#'
#' @description Gaussian : prediction, binomial/multinomial:probability, poisson: expected value
#' @param ypred predicted link function values
#' @param y observations
#' @param family model family
#' @param model0 trained model
deviances <- function(ypred, y, family, model0, weights=NULL){
    if(family == "gaussian"){
        errors = (ypred-y)^2
    }else if (family == "binomial"){
        y0 =model0$glmnet.fit$classnames
        y = ifelse(y == y0[1], 0, 1)
        ypred = ifelse(ypred < 10^(-4), 10^(-4), ypred)
        ypred = ifelse(ypred > 1-10^(-4), 1-10^(-4),  ypred)
        errors = -2*y*log(ypred)-2*(1-y)*log(1-ypred)
    }else if (family == "multinomial"){
        y0 = model0$glmnet.fit$classnames
        YY = matrix(0, ncol = ncol(ypred), nrow(ypred))
        for(i in 1:nrow(YY)){
            YY[i, which(y0 == y[i])] = 1
        }
        ypred = ifelse(ypred < 10^(-4), 10^(-4), ypred)
        ypred = ifelse(ypred > 1-10^(-4), 1-10^(-4),ypred)
        errors = -2*YY*log(ypred)-2*(1-YY)*log(1- ypred)
        errors = apply(errors,1,sum)
    }else if (family == "poisson"){
        errors = 2*(y*log(ypred)-(y-ypred))
    }else{
        if (family == 'cox'){
            errors = logpl1(pred = ypred, surv.time = y[,1], surv.event=y[, 2]) 
        }else {
            stop( paste("model family", family, "not supported!") )
        }
    }
    
    
    if(is.null(weights)  ){ 
     return(errors)
    }else{
        stopifnot(is.vector(weights))
        stopifnot(length(weights) == length(errors))
        weights= length(weights) * weights/sum(weights)
        weights = as.double(weights)
        errors = errors * weights 
        
    }
}


logpl1 <- function(pred, surv.time, surv.event) {	
    
    # from https://github.com/bhklab/survcomp/blob/master/R/logpl.R
    
    n <- length(pred)
    r <- rank(surv.time)
    ita <- pred
    epita <- exp(ita)
    d <- rep(0, n)
    dono <- rep(0, n)
    # for(i in 1:n) {
    #     d[i] <- sum(surv.event[r == r[i]])
    #     dono[i] <- sum(epita[r >= r[i]])
    # }
    
    dono <- sapply(r, function(r_i) sum(epita[r >= r_i]))
    
    if (F) { # wei disabled this: not relevant here.
    risk <- d/dono
    risk1 <- d/dono^{	2}
    culrisk1 <- culrisk <- rep(0, n)
    for(i in 1:n) {
        culrisk[i] <- sum(unique(risk[r <= r[i]]))
        culrisk1[i] <- sum(unique(risk1[r <= r[i]]))
    }
    

    lik <- sum((ita - log(dono)) * surv.event)
    res <- c(lik, sum(surv.event))
    names(res) <- c("logpl", "event")
    return(res)
    }
    
    return( (ita - log(dono)) * surv.event * -2) # need the likelihood for indivual patients. multiply -2 so that the smaller the value the better
}
