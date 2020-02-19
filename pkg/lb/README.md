Conic Log-Binomial Regression
=============================

This **README** provides some additional details on the **lb** package.
The **lb** package was compiled to to further the reproducibility of the
`Conic Log-Binomial Regression` paper.

Simulations A
-------------

The simulation settings A can be reproducted by using the function
`simdata_blizhosm2006`.

    set.seed(0)
    library(lb)

    nobs <- 500L
    setting <- 1L
    d <- simdata_blizhosm2006(nobs, setting = setting)

    s <- lb_convex(d$x, d$y)
    betahat_1 <- setNames(as.double(s), sprintf("beta%i", seq_along(s) - 1L))
    betahat_1

    ##      beta0      beta1 
    ## -2.1627663  0.3619097

    start <- lb:::lb_start(d$y, d$x)
    s <- lb_coptim(d$x, d$y, start = start)
    betahat_2 <- setNames(s$par, sprintf("beta%i", seq_along(s$par) - 1L))
    betahat_2

    ##     beta0     beta1 
    ## -2.163149  0.361980

    s <- lb_auglag(d$x, d$y, start = start)
    betahat_3 <- setNames(s$par, sprintf("beta%i", seq_along(s$par) - 1L))
    betahat_3

    ##      beta0      beta1 
    ## -2.1627662  0.3619097

Simulations B
-------------

The simulation settings B can be reproducted by using the function
`simdata_savuliuyasui2009`.

    nobs <- 500L
    setting <- 1L
    d <- simdata_savuliuyasui2009(nobs, setting = setting)

    s <- lb_convex(d$x, d$y)
    betahat_1 <- setNames(as.double(s), sprintf("beta%i", seq_along(s) - 1L))
    betahat_1

    ##      beta0      beta1      beta2      beta3      beta4      beta5 
    ## -1.8977668  1.1216188 -1.3453598 -0.6935094 -0.7827265 -0.9837980

    start <- lb:::lb_start(d$y, d$x)
    s <- lb_coptim(d$x, d$y, start = start)
    betahat_2 <- setNames(s$par, sprintf("beta%i", seq_along(s$par) - 1L))
    betahat_2

    ##      beta0      beta1      beta2      beta3      beta4      beta5 
    ## -1.8971193  1.1214964 -1.3441929 -0.6942659 -0.7824017 -0.9831321

    s <- lb_auglag(d$x, d$y, start = start)
    betahat_3 <- setNames(s$par, sprintf("beta%i", seq_along(s$par) - 1L))
    betahat_3

    ##      beta0      beta1      beta2      beta3      beta4      beta5 
    ## -1.8977661  1.1216109 -1.3453565 -0.6935069 -0.7827310 -0.9838071

Simulations C
-------------

The simulation settings C can be reproducted by using the function
`sim_data_marschner_plus_minus`.

    nobs <- 500L
    ncov <- 10L
    d <- sim_data_marschner_plus_minus(nobs, ncov)

    s <- lb_convex(d$x, d$y)
    betahat_1 <- setNames(as.double(s), sprintf("beta%i", seq_along(s) - 1L))
    betahat_1

    ##       beta0       beta1       beta2       beta3       beta4       beta5 
    ## -0.54084337 -0.26363609 -0.26363609 -0.13181805 -0.13181805 -0.13181805 
    ##       beta6       beta7       beta8       beta9      beta10 
    ##  0.26363610  0.07725303  0.18638306  0.13181805  0.14538913

    start <- lb:::lb_start(d$y, d$x)
    s <- lb_coptim(d$x, d$y, start = start)
    betahat_2 <- setNames(s$par, sprintf("beta%i", seq_along(s$par) - 1L))
    betahat_2

    ##      beta0      beta1      beta2      beta3      beta4      beta5      beta6 
    ## -0.5409901 -0.2633551 -0.2633552 -0.1316535 -0.1317025 -0.1317017  0.2633552 
    ##      beta7      beta8      beta9     beta10 
    ##  0.0785560  0.1847993  0.1316536  0.1459811

    s <- lb_auglag(d$x, d$y, start = start)
    betahat_3 <- setNames(s$par, sprintf("beta%i", seq_along(s$par) - 1L))
    betahat_3

    ##       beta0       beta1       beta2       beta3       beta4       beta5 
    ## -0.54088185 -0.26361292 -0.26361291 -0.13180637 -0.13180644 -0.13180645 
    ##       beta6       beta7       beta8       beta9      beta10 
    ##  0.26361287  0.07722508  0.18638777  0.13180639  0.14546249
