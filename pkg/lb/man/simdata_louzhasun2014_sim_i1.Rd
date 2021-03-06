% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{simdata_louzhasun2014_sim_i1}
\alias{simdata_louzhasun2014_sim_i1}
\title{Simulate Data based on Luo amd Zhang and Sun 2014}
\usage{
simdata_louzhasun2014_sim_i1(nobs, setting, rtype = c("list",
  "data.frame"))
}
\arguments{
\item{nobs}{number of observations.}

\item{setting}{an integer giving the simulation setting. Allowed values are \eqn{\{1, 2, 3, 4\}}.}

\item{rtype}{a character giving the return type. Allowed values are \code{"list"} and \code{"data.frame"}.}
}
\value{
depending on rtype either a \code{list} or a \code{data.frame}.
}
\description{
Simulate Data based on Luo amd Zhang and Sun 2014
}
\examples{
d <- simdata_louzhasun2014_sim_i1(nobs = 200, 1L)
}
\references{
Ji Luo and Jiajia Zhang and Han Sun. "Estimation of relative risk using a log-binomial model with constraints." Computational Statistics 29.5 (2014): 981-1003.
}
