% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nextdoor.R
\name{nextdoor}
\alias{nextdoor}
\title{Perform the model selection, unbiased error estimation and the nextdoor test(p-value)}
\usage{
nextdoor(errors0, errors, S, nams = NULL, K = 100, B = 1000,
  alpha = 0.1, epsilon = 0.05^2, epsilon2 = 0.05^2, Bindex = NULL,
  pv = TRUE, rescale = TRUE, selectionType = 0, one_sds = rep(0,
  ncol(errors0)), trace = T)
}
\arguments{
\item{errors0}{the n by m errors of the original model sequence.}

\item{errors}{a list of n by m errors of the proximal model sequences.}

\item{S}{a vector corresponding of indexes of the proximal sequences.}

\item{nams}{a length p vector containing feature names. By default, nams = NULL, the features are named according to their column indexes in x.}

\item{K}{number of repetitions estimating the de-biased error}

\item{B}{number of bootstrap repetitions}

\item{alpha}{added error with covariance structure*alpha, by default, alpha = .1}

\item{epsilon}{added error with covariance structure being identity times min(covariance diagonal) times epsilon,
by default, epsilon = 0.05^2}

\item{epsilon2}{added error in the Bootstrap step being min(covariance diagonal) times epsilon2,
by default, epsilon2 = 0.05^2}

\item{Bindex}{n by B index matrix for bootstrap. if Bindex == NULL,
use the default Bootstrap, otherwise, use the provide matrix.}

\item{pv}{if pv == True, estimate the p-values}

\item{rescale}{if rescale == True, perform the mean-rescaled Bootstrap}

\item{selectionType}{if selectionType == 0, pick the model with the smallest randomized error
if selectionType == 1, use the 1se rule}

\item{one_sds}{if the selectionType is 1, the we choose the model with smallest index
such that model error(randomized) <= one_sds[i] + min error(randomized)}

\item{trace}{if trace == True, print the p-value process}
}
\value{
debiased_errors0: de-biased estimate of the model error for the original procedure

debiased_errors: de-biased estimate of the model error for the process excluding a specific feature

worsen: estimated increase in prediction error

pv_value: p values for proximity analysis

result_table: organized result table.
}
\description{
Perform the model selection, unbiased error estimation and the nextdoor test(p-value)
}
\examples{
data(prostateCancerData);library(ranger); library(MASS)
data_train = prostateCancerData$train;nams=prostateCancerData$train$names
data = data.frame(cbind(data_train$feature, data_train$response)); colnames(data) = c(nams, "response")
###get the predictions using the original models and using the nextdoor models with different hyperparameters.###########
###here, for illustration purpose and simplicity, we used the errors from random forest without cross-validation.###########
set.seed(48)
alphas = c(1:10)/10; 
errors0 = array(NA, dim = c(length(data$response),length(alphas))); errors=list()
for(i in 1:length(alphas)){
    model0 = ranger(response~., data = data, alpha = alphas[i])
    errors0[,i] = (data$response-model0$predictions)^2
    for(j in 1:length(nams)){
       if(is.null(errors[j][[1]])){errors[[j]] = array(0, dim = c(length(data$response),length(alphas)))}
       data1 = data[,-j];model0 = ranger(response~., data = data1, alpha = alphas[i]) 
       errors[[j]][,i] = (data$response-model0$predictions)^2
    }
}
res = nextdoor(errors0 = errors0, errors = errors, S =c(1:length(nams)), nams=nams, B = 1000, alpha = 0.1, pv = TRUE, rescale = TRUE, selectionType = 0,trace = TRUE)
print(res, digits =3)
}
