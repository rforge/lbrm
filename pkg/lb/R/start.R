

#' @title Find a starting value
#' @param y vector of observations of length \code{n}.
#' @param x design matrix of dimension \code{n * p}.
#' @param method a character string giving the method used to select the starting values.
#'      The use of \code{"simple"} or \code{"lp"} is recommended.
#' @param delta a numeric additionaly substracted from the intercept.
#' @return a vector giving the starting values.
#' @export
lb_start <- function(y, x, method = c("lp", "simple", "poisson", "unconstrained"), delta = 1) {
    method <- match.arg(method)
    start <- c(-1, double(ncol(x) - 1L))
    if ( method == "lp" ) {
        o <- OP(double(ncol(x)), maximum = TRUE,
                L_constraint(x, leq(nrow(x)), rep(-1e-3, nrow(x))))
        for ( ul_bound in  c(3, 10, 100, Inf)) {
            bounds(o) <- V_bound(ld = -ul_bound, ud = ul_bound, nobj = ncol(x))
            s <- try(ROI_solve(o, "ecos"), silent = TRUE)
            if ( inherits(s, "OP_solution") ) {
                if ( solution(s, "status_code") == 0L ) {
                    start <- solution(s, force = TRUE)
                    break
                }
            }
        }
    } else if ( method == "simple" ) {
        ## used in 
        ## title: Estimation of relative risk using a log-binomial model with constraints
        ## journal: Computational Statistics
        ## author: Luo, Ji and Zhang, Jiajia and Sun, Han
        ## year: 2014
        start <- c(-1, double(ncol(x) - 1L))
    } else if ( method == "poisson" ) {
        ## used in the in lbreg package
        ## title:   Some results for maximum likelihood estimation of adjusted relative risks
        ## journal: Communications in Statistics - Theory and Methods
        ## author:  Bernardo Borba de Andrade and Joanlise Marco de Leon Andrade
        ## year:    2018
        b0 <- coef(glm.fit(x, y, family = poisson(link = "log")))
        mX <- as.matrix(-x[,-1])
        b0[1] <-  min( mX %*% b0[-1] ) - delta
        start <- b0
    } else if ( method == "unconstrained" ) {
        ## The idea is quite simple.
        ## To obtain a starting value we solve the unconstrained problem
        ## a) all the constraints are fullfilled => we are done!
        ## b) we choose a starting value which fullfills all the constraints
        ##    by choosing beta0 small enougth.
        loglike <- function(beta) {
            xb <- drop(x %*% beta)
            -sum(y * xb + (1 - y) * log(1 - exp(xb)))
        }
        gradient <- function(beta) {
            exb <- exp(drop(x %*% beta))
            -drop(crossprod(x, (y - exb) / (1 - exb)))
        }
        start <- c(-1, double(ncol(x) - 1))
        cntrl <- list(fnscale = -1, maxit = 100)
        b0 <- optim(start, loglike, gradient, method = "BFGS", control = cntrl)$par
        mX <- as.matrix(-x[,-1])
        b0[1] <-  min( mX %*% b0[-1] ) - delta
        start <- b0
    }
    return(start)
}

