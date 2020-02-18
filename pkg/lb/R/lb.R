## Imports
#' @importFrom stats rnorm qnorm dbinom rbinom runif binomial coef constrOptim optim poisson setNames terms glm.fit
#' @importFrom utils head modifyList
#' @importFrom slam simple_triplet_zero_matrix simple_triplet_diag_matrix simple_triplet_matrix
#' @importFrom alabama auglag
#' @importFrom R.methodsS3 throw
#' @importFrom R.oo getMessage
#' @importFrom R.utils TimeoutException Arguments
#' @import ROI
#' @import ROI.plugin.ecos
#' @import ROI.plugin.scs
#' @import ROI.plugin.lpsolve


stdm <- simple_triplet_diag_matrix
stzm <- simple_triplet_zero_matrix

lb_convex_dense <- function(x, y, tol = 1e-8, ceps = 1e-7, control = list(), 
    dry_run = FALSE, solver = "ecos") {

    y_is_0 <- y == 0L
    n_y_is_0 <- sum(y_is_0)
    o <- OP(c(y %*% x, double(n_y_is_0), rep(1, n_y_is_0)), maximum = TRUE)
    L1 <- cbind(x, matrix(0, nrow(x), 2 * n_y_is_0))
    log1exp <- function(xi, j, n_y_is_0) {
        M <- matrix(0, nrow = 6, ncol = length(xi) + 2 * n_y_is_0)
        M[1, seq_along(xi)] <- -xi
        M[3, length(xi) + j] <- -1
        M[4, length(xi) + n_y_is_0 + j] <- -1
        M[6, length(xi) + j] <- 1
        M
    }
    L2 <- mapply(log1exp, split(x[y_is_0,], seq_len(n_y_is_0)),
        seq_len(n_y_is_0), MoreArgs = list(n_y_is_0 = n_y_is_0),
        SIMPLIFY = FALSE)
    rhs <- c(c(0, 1, 0), c(0, 1, 1))
    rhs <- c(rep(-ceps, nrow(x)), rep(rhs, n_y_is_0))
    cones <- c(K_lin(nrow(x)), K_expp(2 * n_y_is_0))
    L <- do.call(rbind, c(list(L1), L2))
    constraints(o) <- C_constraint(L, cones, rhs)
    bounds(o) <- V_bound(ld = -Inf, nobj = length(objective(o)))
    if ( dry_run ) return(o)
    control$tol <- tol
    so <- ROI_solve(o, solver = solver, control = control)
    rval <- head(solution(so, force = TRUE), NCOL(x))
    if ( solution(so, "status_code") != 0L ) {
        warning("Solution: ", paste(solution(so, "status")$msg$message, collapse = " "))
    }
    attributes(rval)$solution <- so
    rval
}


lb_convex_sparse <- function(x, y, tol = 1e-8, ceps = 1e-7, control = list(), 
    dry_run = FALSE, solver = "ecos") {

    y_is_0 <- y == 0L
    n_y_is_0 <- sum(y_is_0)
    o <- OP(c(y %*% x, double(n_y_is_0), rep(1, n_y_is_0)), maximum = TRUE)
    L1 <- cbind(x, stzm(nrow(x), 2 * n_y_is_0))
    log1exp <- function(xi, j, n_y_is_0) {
        k <- length(xi)        
        si <- c(rep.int(1, k),     3,                4,     6)
        sj <- c(seq_along(xi), k + j, k + n_y_is_0 + j, k + j)
        sv <- c(          -xi,    -1,               -1,     1)
        simple_triplet_matrix(si, sj, sv, nrow = 6L, ncol = (k + 2 * n_y_is_0))
    }
    L2 <- mapply(log1exp, split(x[y_is_0,], seq_len(n_y_is_0)),
        seq_len(n_y_is_0), MoreArgs = list(n_y_is_0 = n_y_is_0),
        SIMPLIFY = FALSE)
    rhs <- c(c(0, 1, 0), c(0, 1, 1))
    rhs <- c(rep(-ceps, nrow(x)), rep(rhs, n_y_is_0))
    cones <- c(K_lin(nrow(x)), K_expp(2 * n_y_is_0))
    L <- do.call(rbind, c(list(L1), L2))
    constraints(o) <- C_constraint(L, cones, rhs)
    bounds(o) <- V_bound(ld = -Inf, nobj = length(objective(o)))
    if ( dry_run ) return(o)
    control$tol <- tol
    so <- ROI_solve(o, solver = solver, control = control)
    rval <- head(solution(so, force = TRUE), NCOL(x))
    if ( solution(so, "status_code") != 0L ) {
        warning("Solution: ", paste(solution(so, "status")$msg$message, collapse = " "))
    }
    attributes(rval)$solution <- so
    rval
}


order_matrix <- function(x) {
    do.call(order, as.list(as.data.frame(x)))
}

group_matrix <- function(x) {
    i <- order_matrix(x)
    x <- x[i,]
    b <- duplicated(x)
    group_index <- cumsum(!b)
    x <- x[!b,]
    attributes(x)$weights <- tabulate(group_index)
    x
}

lb_convex_weighted <- function(x, y, tol = 1e-8, ceps = 1e-7, control = list(), 
    dry_run = FALSE, solver = "ecos") {

    y_is_0 <- y == 0L
    x0 <- x[y_is_0,]
    x0 <- group_matrix(x0)
    x1 <- unique(x[!y_is_0,])
    n <- nrow(x0) + nrow(x1)
    o <- OP(c(y %*% x, double(nrow(x0)), attr(x0, "weights")), maximum = TRUE)
    L1 <- cbind(rbind(x0, x1), matrix(0, n, 2 * nrow(x0)))
    log1exp <- function(xi, j, n_y_is_0) {
        M <- matrix(0, nrow = 6, ncol = length(xi) + 2 * n_y_is_0)
        M[1, seq_along(xi)] <- -xi
        M[3, length(xi) + j] <- -1
        M[4, length(xi) + n_y_is_0 + j] <- -1
        M[6, length(xi) + j] <- 1
        M
    }
    L2 <- mapply(log1exp, split(x0, seq_len(nrow(x0))),
        seq_len(nrow(x0)), MoreArgs = list(n_y_is_0 = nrow(x0)),
        SIMPLIFY = FALSE)
    rhs <- c(c(0, 1, 0), c(0, 1, 1))
    rhs <- c(rep(-ceps, n), rep(rhs, nrow(x0)))
    cones <- c(K_lin(n), K_expp(2 * nrow(x0)))
    L <- do.call(rbind, c(list(L1), L2))
    constraints(o) <- C_constraint(L, cones, rhs)
    bounds(o) <- V_bound(ld = -Inf, nobj = length(objective(o)))
    if ( dry_run ) return(o)
    control$tol <- tol
    so <- ROI_solve(o, solver = solver, control = control)
    rval <- head(solution(so, force = TRUE), NCOL(x))
    if ( solution(so, "status_code") != 0L ) {
        warning("Solution: ", paste(solution(so, "status")$msg$message, collapse = " "))
    }
    attributes(rval)$solution <- so
    return(rval)
}


#' @title Fitting Log-Binomial Models by Convex Optimization
#' @param x design matrix of dimension \code{n * p}.
#' @param y vector of observations of length \code{n}.
#' @param method a character string giving the method used to construct the optimization problem,
#'     possible
#' @param tol tolerance for the optimizer (default is \code{1e-8}).
#' @param ceps epsilon subtracted from the right hand side (\eqn{X \beta \leq 0 - ceps}).
#'  Since the inequality \eqn{X \beta \leq 0} is only fullfilled to a certain tolerance
#'  we have to substact epsilon to ensure \eqn{X \beta \leq 0}.
#' @param control a list containing additional arguments passed to \code{ROI_solve}.
#' @param dry_run a logical controling if the problem should be solved 
#'     (default \code{dry\_run = FALSE}) or the optimization problem should be returned.
#' @param solver a character string selecting the solver to be used (default is \code{"ecos"}).
#' @return a vector giving the estimated coefficients.
#' @examples
#' d <- simdata_blizhosm2006(500, 1)
#' lb_convex(d$x, d$y)
#' @export
lb_convex <- function(x, y, method = c("dense", "sparse", "weighted"), 
    tol = 1e-8, ceps = 1e-7, control = list(), dry_run = FALSE, solver = "ecos") {

    method <- match.arg(method)

    if ( isTRUE(method == "dense") ) {
        s <- lb_convex_dense(x = x, y = y, tol = tol, ceps = ceps, control = control, 
                             dry_run = dry_run, solver = solver)
    } else if ( isTRUE(method == "sparse") ) {
        s <- lb_convex_sparse(x = x, y = y, tol = tol, ceps = ceps, control = control, 
                              dry_run = dry_run, solver = solver)
    } else if ( isTRUE(method == "weighted") ) {
        s <- lb_convex_weighted(x = x, y = y, tol = tol, ceps = ceps, control = control, 
                                dry_run = dry_run, solver = solver)
    } else {
        stop("not implementated!")
    }
    return(s)
}


neg_log_likelihood_default <- function(x, y, beta) {
    eta <- drop(tcrossprod(beta, x))
    -sum(dbinom(x = y, size = 1L, prob = exp(eta), log = TRUE))
}

elog <- function(x) ifelse(x > 0, log(x), -Inf)

neg_log_likelihood_auglag <- function(x, y, beta) {
    eta <- drop(tcrossprod(beta, x))
    b <- (y == 1)
    -sum(eta[b]) - sum(elog(1 - exp(eta[!b])))
}

gradient_default <- function(x, y, beta) {
    mu <- exp(drop(tcrossprod(beta, x)))
    -drop(crossprod(x, (y - mu) / (1 - mu)))
}

gradient_auglag <- function(x, y, beta) {
    eta <- drop(tcrossprod(beta, x))
    mu <- pmin(exp(eta), 1 - .Machine$double.neg.eps)
    -drop(crossprod(x, (y - mu) / (1 - mu)))
}


#' @title Fitting Log-Binomial Models with \code{constrOptim}
#' @param x design matrix of dimension \code{n * p}.
#' @param y vector of observations of length \code{n}.
#' @param start numeric (vector) starting value (of length p): must be in the
#'              feasible region.
#' @param mu see \code{constrOptim}.
#' @param control passed to ``optim''.
#' @param method passed to ``optim''.
#' @param outer_iterations see \code{constrOptim}.
#' @param outer_eps see \code{constrOptim}.
#' @param ... see \code{constrOptim}.
#' @return a vector giving the estimated coefficients.
#' @export
lb_coptim <- function(x, y, start, mu = 1e-04, control = list(maxit = 500L), 
    method = "BFGS", outer_iterations = 10000, outer_eps = 1e-08, ...) {

    stopifnot(length(start) == ncol(x))

    nloglike <- function(beta) neg_log_likelihood_default(x, y, beta)
    gradient <- function(beta) gradient_default(x, y, beta)

    grad <- if (method %in% c("BFGS", "CG", "L-BFGS-B")) gradient else NULL
    s <- constrOptim(start, f = nloglike, grad = grad, ui = -x, ci = 0, mu = mu,
        control = control, method = method, outer.iterations = outer_iterations, 
        outer.eps = outer_eps, ..., hessian = FALSE)
    return(s)
}


#' @title Fitting Log-Binomial Models with \code{auglag} from the \pkg{alabama} Package
#' @param x design matrix of dimension \code{n * p}.
#' @param y vector of observations of length \code{n}.
#' @param start numeric (vector) starting value (of length p).
#' @param tol a numeric giving the tolerance for convergence of the outer iterations 
#'     (see \code{eps} in \code{control.outer} the from the \code{auglag} function.
#' @param ceps epsilon subtracted from the right hand side (\eqn{X \beta \leq 0 - ceps}).
#'     Since the inequality \eqn{X \beta \leq 0} is only fullfilled to a certain tolerance
#'     we have to substact epsilon to ensure \eqn{X \beta \leq 0}.
#' @param control passed to \code{control.outer}.
#' @param control.optim passed to \code{control.optim}.
#' @param implementation a character string choosing which implementation of the 
#'     log-likelihood and gradient is used.
#' @return the optimization result.
#' @export
lb_auglag <- function(x, y, start, tol = 1e-7, ceps = 1e-7, control = list(), 
    control.optim = list(), implementation = c("improved", "naive")) {

    stopifnot(length(start) == ncol(x))
    implementation <-  match.arg(implementation)

    if ( implementation == "naive" ) {
        nloglike <- function(beta) neg_log_likelihood_default(x, y, beta)
        gradient <- function(beta) gradient_default(x, y, beta)
    } else {
        nloglike <- function(beta) neg_log_likelihood_auglag(x, y, beta)
        gradient <- function(beta) gradient_auglag(x, y, beta)
    }
    
    constraint <- function(beta) -drop(x %*% beta) - ceps
    hin.jac <- function(beta) -x

    control$eps <- tol
    control$trace <- FALSE

    s <- auglag(par = start, fn = nloglike, gr = gradient, 
      hin = constraint, hin.jac = hin.jac, 
      control.outer = control, control.optim = control.optim)
    return(s)
}


set_default_options <- function(x) {
    default_options <- list(mu_init = 0.1, print_level = 0, max_iter = 5000, tol = 1e-8,
        jac_c_constant = "yes", jac_d_constant = "yes", linear_solver = "ma97")
    n <- names(default_options)[!names(default_options) %in% names(x)]
    x[n] <- default_options[n]
    x
}


#' @title Fitting Log-Binomial Models with \code{ipoptr}
#' @param x design matrix of dimension \code{n * p}.
#' @param y vector of observations of length \code{n}.
#' @param start numeric (vector) starting value (of length p).
#' @param tol a numeric giving the tolerance for convergence of the outer iterations 
#'     (see \code{eps} in \code{control.outer} the from the \code{auglag} function.
#' @param ceps epsilon subtracted from the right hand side (\eqn{X \beta \leq 0 - ceps}).
#'     Since the inequality \eqn{X \beta \leq 0} is only fullfilled to a certain tolerance
#'     we have to substact epsilon to ensure \eqn{X \beta \leq 0}.
#' @param control passed to \code{control.outer}.
#' @return the optimization result.
#' @export
lb_ipopt <- function(x, y, start, tol = 1e-8, ceps = 1e-7, control = list()) {
    
    stopifnot(length(start) == ncol(x))

    suppressWarnings(is_ipopts_installed <- requireNamespace("ipoptr", quietly = TRUE))
    if (!is_ipopts_installed) stop("ipoptr is not installed!")

    nloglike <- function(beta) neg_log_likelihood_default(x, y, beta)
    gradient <- function(beta) gradient_default(x, y, beta)

    hessi <- function(x, y, beta) {
        X0 <- x[y == 0L, ]
        p0 <- exp(drop(X0 %*% beta))
        -t(X0) %*% diag(p0 / (1 - p0)^2) %*% X0
    }

    hessian <- function(beta, obj_factor, hessian_lambda) {
        H <- hessi(x, y, beta)
        obj_factor * t(H)[upper.tri(H, TRUE)]
    }

    constraint <- function(beta) {
        drop(x %*% beta)
    }

    lil_matrix <- function(x) {
          y <- unname(x)
          b <- y != 0
          list(which(b, useNames = FALSE), y[b])
    }
    lil <- apply(x, 1, lil_matrix)
    eval_jac_g_structure <- lapply(lil, "[[", 1L)

    eval_jac_g <- function(beta) {
      unlist(lapply(lil, "[[", 2L), FALSE, FALSE)
    }

    control$tol <- tol
    control <- set_default_options(control)

    eval_h_structure <- lapply(seq_len(ncol(x)), seq_len)

    s <- ipoptr::ipoptr(x0 = start, eval_f = nloglike, eval_grad_f = gradient, 
        eval_h = hessian, eval_h_structure = eval_h_structure,
        constraint_lb = rep.int(-1e16, nrow(x)), 
        constraint_ub = rep.int(-ceps, nrow(x)),
        eval_g = constraint, eval_jac_g = eval_jac_g,
        eval_jac_g_structure = eval_jac_g_structure, 
        lb = rep.int(-1e16, ncol(x)), ub = rep.int(1e16, ncol(x)),
        opts = control)
    return(s)
}



