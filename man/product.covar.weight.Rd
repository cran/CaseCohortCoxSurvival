
\name{product.covar.weight}
\alias{product.covar.weight}
\title{product.covar.weight}

\description{
Computes the product of joint design weights and joint sampling indicators 
covariances, needed for the phase-two component of the variance (with design or 
calibrated weights).
}

\usage{
product.covar.weight(casecohort, stratified = NULL)
}

\arguments{
 \item{casecohort}{if \code{stratified = TRUE}, data frame with \code{status} (case 
status), \code{W} (the \eqn{J} strata), \code{strata.m} (vector of length \eqn{J} with the
numbers of sampled individuals in the strata) and \code{strata.n} (vector of length 
\eqn{J} with the strata sizes), for each individual in the stratified case-cohort 
data. If \code{stratified = FALSE}, data frame with \code{status} (case status), 
\code{m} (number of sampled individuals) and \code{n} (cohort size), for each 
individual in the un-stratified case-cohort data.}

  \item{stratified}{was the sampling of the case-cohort stratified on \code{W}? 
Default is \code{FALSE}.}
}

\value{

 \code{product.covar.weight}: matrix with the products of joint design weights and
joint sampling indicator covariances, for the non-cases in the case-cohort.
}

\details{

 \code{product.covar.weight} creates the matrix with the products of joint design 
weights and joint sampling indicator covariances, for the non-cases in the case
cohort. In other words, it has as many rows and columns as non-cases in the case
cohort, and contains the \eqn{w_{i,k,j} \sigma_{i,k,j}}, with 

 \eqn{w_{i,k,j} = \frac{n^{(j)}(n^{(j)} -1)}{m^{(j)}(m^{(j)} -1)}} if individuals
  \eqn{i} and \eqn{k} in stratum \eqn{j} are both non-cases, and 
 \eqn{w_{i,k,j} = \left( \frac{n^{(j)}}{m^{(j)}} \right)^2} otherwise, 
 \eqn{i \neq k \in \lbrace 1, \dots, n^{(j)} \rbrace},
 \eqn{j \in \lbrace 1, \dots, J \rbrace}.


 \eqn{w_{i,i,j} = \frac{n^{(j)}}{m^{(j)}}} if individuals \eqn{i} in stratum \eqn{j} 
  is a non-case, \eqn{i \in \lbrace 1, \dots, n^{(j)} \rbrace}, 
   \eqn{j \in \lbrace 1, \dots, J \rbrace}.

 \eqn{\sigma_{i,k,j} = \frac{m^{(j)}(m^{(j)} -1)}{n^{(j)}(n^{(j)} -1)} - 
    \left( \frac{m^{(j)}}{n^{(j)}} \right)^2} if individuals \eqn{i} and
   \eqn{k} in stratum \eqn{j} are both non-cases, 
  \eqn{i \neq k \in \lbrace 1, \dots, n^{(j)} \rbrace},
   \eqn{j \in \lbrace 1, \dots, J \rbrace}.

 \eqn{\sigma_{i,i,j} = \frac{m^{(j)}}{n^{(j)}} - \left(1 - \frac{m^{(j)}}{n^{(j)}} \right)}
   if individuals \eqn{i} in stratum \eqn{j} is a non-case, 
  \eqn{i \in \lbrace 1, \dots, n^{(j)} \rbrace}, 
  \eqn{j \in \lbrace 1, \dots, J \rbrace}.

See Section 3.3 in Etievant and Gail (2023).
}

\references{
Etievant, L., Gail, M.H. (2023). Cox model inference for relative hazard and 
pure risk from stratified weight-calibrated case-cohort data. Submitted.
}

\seealso{

\code{\link{variance}}, that uses \code{product.covar.weight} to compute the variance 
estimate that follows the complete variance decomposition (superpopulation and 
phase-two variance components).
}

\examples{
data(dataexample, package="CaseCohortCoxSurvival")

casecohort  <- dataexample$casecohort # a simulated stratified case-cohort

prod.covar.weight <- product.covar.weight(casecohort, stratified = TRUE)

nrow(prod.covar.weight)
ncol(prod.covar.weight)
sum(casecohort$status == 0) # number of non-cases in the case-cohort
}
