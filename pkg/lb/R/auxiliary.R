
timeit <- function(expr, substitute = TRUE, envir = parent.frame(), timeout, 
    cpu = timeout, elapsed = timeout) {
    time_stop <- time_start <- NA_real_

    if (substitute) expr <- substitute(expr)
    if (!is.environment(envir)) throw("Argument 'envir' is not a list: ", class(envir)[1L])

    cpu <- Arguments$getNumeric(cpu, range = c(0, Inf))
    elapsed <- Arguments$getNumeric(elapsed, range = c(0, Inf))
    
    setTimeLimit(cpu = cpu, elapsed = elapsed, transient = TRUE)
    on.exit({
        setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    })
    ret <- tryCatch({
        time_start <- Sys.time()
        ret <- eval(expr, envir = envir)
        time_stop <- Sys.time()
        ret
    }, error = function(ex) {
        msg <- ex$message
        pattern <- gettext("reached elapsed time limit", "reached CPU time limit", domain = "R")
        pattern <- paste(pattern, collapse = "|")
        if (regexpr(pattern, msg) != -1L) {
            ex <- TimeoutException(msg, cpu = cpu, elapsed = elapsed)
            warning(getMessage(ex))
            structure(list(), class = c("timeout", "list"))
        } else {
            throw(ex)
        }
    })
    attr(ret, "time") <- list(start = time_start, stop = time_stop, elapsed = if(is.na(time_stop)) NA_real_ else time_stop - time_start)
    ret
}


lb_check_constraints <- function(x, beta) {
	constraints_fullfilled <- drop(x %*% beta) <= 0
    
    if (all(constraints_fullfilled)) return(TRUE)
    
    rval <- FALSE
    tab <- table(constraints_fullfilled)
    ptab <- prop.table(tab)
    n_constr_violations <- tab["FALSE"]
    percent_constr_violations <- 100 * ptab["FALSE"]
    attributes(rval)$message <- sprintf("%i constraint violations (%0.2f %%)", 
        n_constr_violations, percent_constr_violations)
    attributes(rval)$constr_fullfilled <- tab
    return(rval)
}

