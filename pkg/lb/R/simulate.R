
simlist <- function(...) {
    x <- list(...)
    class(x) <- c("simlist", class(x))
    x
}

as.simlist <- function(x) UseMethod("as.simlist", x)

as.simlist.list <- function(x) {
    class(x) <- unique(c("simlist", class(x)))
    x
}

#' @noRd
#' @export
print.simlist <- function(x, ...) {
    class(x) <- class(x)[class(x) != "simlist"]
    print(x)
}

## we assume that has always a intercept
#' @noRd
#' @export
as.data.frame.simlist <- function(x, ...) {
    cbind(data.frame(y = x$y), x$x[,-1])
}


#' @title Simulate Data based on Blizzard and Hosmer 2006
#' @param nobs number of observations.
#' @param setting an integer giving the simulation setting. Allowed values are 
#'      \eqn{\{1, 2, 3, 4, 5, 6, 7, 8\}}.
#' @param rtype a character giving the return type. Allowed values are \code{"list"} and \code{"data.frame"}.
#' @return depending on rtype either a \code{list} or a \code{data.frame}.
#' @references
#' L. Blizzard and D. W. Hosmer. "Parameter estimation and goodness-of-fit in log binomial regression." Biometrical Journal 48.1 (2006): 5-22.
#' @examples
#' d <- simdata_blizhosm2006(500, 1)
#' @export
simdata_blizhosm2006 <- function(nobs, setting, rtype = c("list", "data.frame")) {
    stopifnot(isTRUE(setting %in% 1:8))
    rtype <- match.arg(rtype)
    k <- as.integer(setting)
    b0 <- rep(c(-2.30259, -1.20397, -0.69315, -0.35667), each = 2)
    b1 <- rep(c(0.38376, 0.56687, 0.65200, 0.70808), each = 2)
    a <- c(6, 4, 2, 1, 1, 0, 0.5, -0.5)
    x <- runif(nobs, -6, a[k])
    yhat <- b0[k] + b1[k] * x
    y <- rbinom(length(x), 1, exp(b0[k] + b1[k] * x))
    rval <- if (rtype == "list") simlist(y = y, x = cbind(intercept = 1, x)) else data.frame(y = y, x = x)
    attributes(rval)$yhat <- yhat
    attributes(rval)$b0 <- b0[k]
    attributes(rval)$b1 <- b1[k]
    attributes(rval)$a <- a[k]
    rval
}

simdata_blizhosm2006_2 <- function(nobs, setting, rtype = c("list", "data.frame")) {
    rtype <- match.arg(rtype)
    stopifnot(setting %in% 1:4)
    params <- list(c(log(0.3), log(1.5), 0.18, 0.2), c(log(0.3), log(1.5), 0.18, 0.5),
                   c(log(0.3), log(2.0), 0.10, 0.2), c(log(0.3), log(2.0), 0.10, 0.5))
    prob <- params[[setting]][4L]
    beta <- params[[setting]][-4L]
    d <- .simdata_louzhasun2014_sim_i1(nobs, prob, beta, rtype)
    attr(d, "beta") <- beta
    d
}


#' @title Simulate Data similar to an example in the \pkg{rms} Manual
#' @param nobs number of observations.
#' @param rtype a character giving the return type. Allowed values are \code{"list"} and \code{"data.frame"}.
#' @return depending on rtype either a \code{list} or a \code{data.frame}.
#' @references
#' Frank E Harrell Jr (2019). rms: Regression Modeling Strategies. R package version 5.1-3. https://CRAN.R-project.org/package=rms
#' \url{https://www.rdocumentation.org/packages/rms/versions/5.1-2/topics/anova.rms}
#' @examples
#' d <- simdata_rms(500)
#' @export
simdata_rms <- function(nobs, rtype = c("list", "data.frame")) {
    rtype <- match.arg(rtype)
    treat <- factor(sample(c('a','b','c'), size = nobs, replace = TRUE))
    num.diseases <- sample(0:4, nobs, replace = TRUE)
    age <- rnorm(nobs, 50L, 10L)
    cholesterol <- rnorm(nobs, 200L, 25L)
    weight <- rnorm(nobs, 150L, 20L)
    sex <- factor(sample(c('female','male'), size = nobs, replace = TRUE))

    # Specify population model for log odds that Y = 1
    L <- ( -1 + 0.1 * (num.diseases - 2) + 0.045 * (age - 70) 
         + (log(cholesterol - 10) - 5.2) - 2 * (treat == 'a') 
         + 0.5 * (treat == 'b') - 0.5 * (treat == 'c') )
    # Simulate binary y to have Prob(y = 1) = 1 / (1 + exp(-L))
    y <- as.double(runif(nobs) < exp(L))

    x <- cbind(intercept = 1, age, sex, weight, 
        logchol = log(cholesterol - 10), num.diseases, 
        treatb = (treat == "b"), treatc = (treat == "c"))

    if ( rtype == "list" ) simlist(y = y, x = x) else cbind(data.frame(y = y), x[,-1])
}


#' @title Simulate Data based on Luo amd Zhang and Sun 2014
#' @param nobs number of observations.
#' @param setting an integer giving the simulation setting. Allowed values are \eqn{\{1, 2, 3, 4\}}.
#' @param rtype a character giving the return type. Allowed values are \code{"list"} and \code{"data.frame"}.
#' @return depending on rtype either a \code{list} or a \code{data.frame}.
#' @references
#' Ji Luo and Jiajia Zhang and Han Sun. "Estimation of relative risk using a log-binomial model with constraints." Computational Statistics 29.5 (2014): 981-1003.
#' @examples
#' d <- simdata_louzhasun2014_sim_i1(nobs = 200, 1L)
#' @export
simdata_louzhasun2014_sim_i1 <- function(nobs, setting, rtype = c("list", "data.frame")) {
    rtype <- match.arg(rtype)
    stopifnot(setting %in% 1:4)
    prob <- 0.2
    betas <- list(c(log(0.3), log(2.0), 0.10), c(log(0.3), log(1.5), 0.18),
                  c(log(0.3), log(2.0), 0.12), c(log(0.3), log(1.5), 0.19))
    beta <- betas[[setting]]
    d <- .simdata_louzhasun2014_sim_i1(nobs, prob, beta, rtype)
    attr(d, "beta") <- beta
    d
}

.simdata_louzhasun2014_sim_i1 <- function(nobs, prob, beta, rtype = c("list", "data.frame")) {
    rtype <- match.arg(rtype)
    stopifnot(is.numeric(prob), length(prob) == 1L, prob >= 0, prob <= 1)
    x0 <- rbinom(nobs, 1, prob)
    n0 <- sum(as.numeric(x0 == 0))
    x1 <- double(nobs)
    x1[x0 == 0] <- runif(n0, -6, 2)
    x1[x0 == 1] <- runif(nobs - n0, -4, 4)
    mu <- exp(cbind(1, x0, x1) %*% beta)
    ## NOTE: Their reported code is inconsistent with their paper!
    stopifnot(all(mu <= 1))
    y <- rbinom(length(mu), 1, mu)
    x <- cbind(intercept = 1, x0 = x0, x1 = x1)
    if ( rtype == "list" ) simlist(y = y, x = x) else cbind(data.frame(y = y), x[,-1])
}

simdata_louzhasun2014_sim_i2 <- function(nobs, rtype = c("list", "data.frame")) {
    rtype <- match.arg(rtype)
    prob <- 0.5
    beta <- c(log(0.25), log(2), 0.11, 0.15)
    x0 <- rbinom(nobs, 1, prob)
    n0 <- sum(as.numeric(x0 == 0))
    x1 <- double(nobs)
    x1[x0 == 0] <- runif(n0, -4, 0)
    x1[x0 == 1] <- runif(nobs - n0, -2, 2)
    x2 <- x1 * rnorm(nobs, 1, 0.1)
    mu <- exp(cbind(1, x0, x1, x2) %*% beta)
    stopifnot(all(mu <= 1))
    y <- rbinom(length(mu), 1, mu)
    x <- cbind(intercept = 1, x0 = x0, x1 = x1)
    if ( rtype == "list" ) simlist(y = y, x = x) else cbind(data.frame(y = y), x[,-1])
}


#' @title Simulate Data based on Marschner 2018
#' @param nobs number of observations.
#' @param ncov number of covariates.
#' @param setting an integer giving the simulation setting. Allowed values are \eqn{\{1, 2\}}.
#' @param rtype a character giving the return type. Allowed values are \code{"list"} and \code{"data.frame"}.
#' @return depending on rtype either a \code{list} or a \code{data.frame}.
#' @examples
#' d <- simdata_marschner_jss_2018(nobs = 2000, ncov = 20, setting = 1)
#' @export
simdata_marschner_jss_2018 <- function(nobs, ncov, setting, rtype = c("list", "data.frame")) {
    stopifnot(setting %in% c(1L, 2L))
    rtype <- match.arg(rtype)
    ## Number of covariates
    # ncov <- c(5, 10, 15)
    ## Relative risk associated with each covariate
    rr <- c(0.8, 1)[setting]
    # Baseline risk
    p0 <- 0.6
    x <- matrix(round(runif(nobs * ncov)), nrow = nobs, ncol = ncov)
    colnames(x) <- sprintf("X%i", seq_len(NCOL(x)))
    fitprob <- p0 * exp(rowSums(log(rr) * x))
    y <- rbinom(nobs, 1, fitprob)
    stopifnot(qr(x[y == 0,])$rank == ncol(x))
    if ( rtype == "list" ) simlist(y = y, x = cbind(1, x)) else cbind(data.frame(y = y), x)
}


#' @title Simulate Data based on Savu, Liu and Yasui 2009
#' @param nobs number of observations.
#' @param setting an integer giving the simulation setting. Allowed values are \eqn{\{1, 2, 3, 4, 5, 6, 7, 8\}}.
#' @param rtype a character giving the return type. Allowed values are \code{"list"} and \code{"data.frame"}.
#' @return depending on rtype either a \code{list} or a \code{data.frame}.
#' @references
#' Savu A, Liu Q, Yasui Y (2010) Estimation of relative risk and prevalence ratio. Statistics in Medicine 29(22):2269-2281, DOI \url{10.1002/sim.3989}
#' @examples
#' d <- simdata_savuliuyasui2009(nobs = 200, 1L)
#' @export
simdata_savuliuyasui2009 <- function(nobs, setting, rtype = c("list", "data.frame")) {
    rtype <- match.arg(rtype)
    stopifnot(setting %in% 1:8)
    x1 <- rbinom(nobs, 1, prob = 0.5)
    x2 <- sample(c(0, 1, 2), nobs, TRUE, prob = c(0.3, 0.3, 0.4))
    x21 <- as.integer(x2 == 1)
    x22 <- as.integer(x2 == 2)
    x3 <- runif(nobs, -1, 2)
    beta <- c(-1, 1, -1, 1, 1)
    x <- cbind(1, x1, x21, x22, x3)
    p_exposure <- binomial(link = "logit")$linkinv(drop(x %*% beta))
    exposure <- rbinom(nobs, 1, p_exposure)
    hx <- x1 + x21 + x22 + x3
    b <- exposure == 1L
    fun_exposure0 <- switch(setting,
        "1" = function(x) binomial(link = "log")$linkinv(-2.1 - x),
        "2" = function(x) binomial(link = "cloglog")$linkinv(-1.9 - x),
        "3" = function(x) binomial(link = "logit")$linkinv(-1.7 - x),
        "4" = function(x) binomial(link = "probit")$linkinv(-1.48 - x),
        "5" = function(x) binomial(link = "log")$linkinv(-1.1 - pmax(x, 0)),
        "6" = function(x) binomial(link = "cloglog")$linkinv(-0.9 - pmax(x, 0)),
        "7" = function(x) binomial(link = "logit")$linkinv(-0.7 - pmax(x, 0)),
        "8" = function(x) binomial(link = "probit")$linkinv(-0.48 - pmax(x, 0)))
    fitprob <- fun_exposure0(hx)
    fitprob[b] <- 3 * fitprob[b]
    y <- rbinom(nobs, 1, prob = pmin(fitprob, 1))
    x <- cbind(intercept = 1, exposure = exposure, x1 = x1, x21 = x21, x22 = x22, x3 = x3)
    stopifnot(qr(x[y == 0,])$rank == ncol(x))
    d <- if ( rtype == "list" ) simlist(y = y, x = x) else cbind(data.frame(y = y), x)
    attr(d, "beta") <- beta
    d
}


#' @title Simulate Data similar to Marschner 2018
#' @description This slightly modified verison allows to use more covariates.
#' @param nobs number of observations.
#' @param ncov number of covariates.
#' @return depending on rtype either a \code{list} or a \code{data.frame}.
#' @examples
#' d <- sim_data_marschner_plus_minus(nobs = 200, ncov = 20L)
#' @export
sim_data_marschner_plus_minus <- function(nobs, ncov) {
    stopifnot(abs(ncov / 2 - as.integer(ncov / 2)) < 1e-12)
    .sim_data_marschner_plus_minus <- function(nobs, ncov, rr1, rr2, p0 = 0.6) {
        x <- cbind(1, matrix(round(runif(nobs * ncov)), nrow = nobs, ncol = ncov))
        colnames(x) <- c("intercept", sprintf("x%i", seq_len(NCOL(x) - 1L)))
        beta <- c(log(p0), log(c(rep.int(rr1, ncov / 2), rep.int(rr2, ncov / 2))))
        fitprob <- exp(x %*% beta)
        y <- rbinom(nobs, 1, pmin(fitprob, 1))
        if ( (qr(x[y == 0,])$rank != ncol(x)) | (qr(x[y == 1,])$rank != ncol(x)) )
            return(NULL)
        list(y = y, x = x, beta = beta, fitprob = fitprob, meta = list(rr = c(rr1, rr2), p0 = p0))
    }
    rr1 <- exp(-1)
    rr2 <- exp( 1)
    d <- NULL
    for (i in seq_len(100)) {
        d <- .sim_data_marschner_plus_minus(nobs = nobs, ncov = ncov, rr1 = rr1, rr2 = rr2, p0 = 0.6)
        if ( !is.null(d) ) break
    }
    d$iter <- i
    d
}

