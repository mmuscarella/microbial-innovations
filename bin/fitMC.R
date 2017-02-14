################################################################################
#
# fitMC Function for Ancestral State Reconstruction
#
# Written by: Mario Muscarella
#
# Last Update: 20170213
#
# Goals:
#       1. 
#
################################################################################

fitMC <- function(x, phy, ip = 0.5){

  # Test
  # sim <- TT.sim(birth = 0.2, a = 0.90, b = 0.90, tips = 100)
  # phy <- sim$tree
  # x <- sim$traits

  # Define Tree Features
  nb.tips <- length(phy$tip.label)
  nb.node <- phy$Nnode

  # Define Output Object
  obj <- list()

  # Traits as Factors
  x <- x[phy$tip.label]
  x <- factor(x)
  nl <- nlevels(x)
  lvls <- levels(x)
  x <- as.integer(x)

  # Define Rate Matrix
  rate <- matrix(NA, nl, nl)
  np <- nl * (nl - 1)
  rate[col(rate) != row(rate)] <- 1:np
  index.matrix <- rate
  tmp <- cbind(1:nl, 1:nl)
  index.matrix[tmp] <- NA
  rate[tmp] <- 0
  rate[rate == 0] <- np + 1

  # Define Likelihood Output
  liks <- matrix(0, nb.tips + nb.node, nl)
  TIPS <- 1:nb.tips
  liks[cbind(TIPS, x)] <- 1

  phy <- reorder(phy, "postorder")
  H <- matrix(0, nl, nl)
  e1 <- phy$edge[, 1]
  e2 <- phy$edge[, 2]
  EL <- phy$edge.length

  # Define Function
  E <- expm::expm
  dev <- function(p){
    if (any(is.nan(p)) || any(is.infinite(p)))
      return(1e+50)
    comp <- numeric(nb.tips + nb.node)
    H[] <- c(p, 0)[rate]
    diag(H) <- -rowSums(H)
    for (i in seq(from = 1, by = 2, length.out = nb.node)){
      j <- i + 1L
      anc <- e1[i]
      des1 <- e2[i]
      des2 <- e2[j]
      v.l <- E(H * EL[i]) %*% liks[des1, ]
      v.r <- E(H * EL[j]) %*% liks[des2, ]
      v <- v.l * v.r
      comp[anc] <- sum(v)
      liks[anc, ] <- v/comp[anc]
    }
    #return(liks[-TIPS,])
    #dev <- -2 * sum(log(comp[-TIPS]))
    dev <- -sum(log(comp[-TIPS]))

    return(dev)
  }

  out <- optim(rep(ip, length.out = np), dev, NULL,
        method = "L-BFGS-B",
        lower = rep(1e-2, np), upper = rep(1e+2, np))

  #obj$loglik <- -out$value/2
  #obj$rates <- out$par

  #out.nlm <- nlm(function(p) dev(p), p = obj$rates, iterlim = 1,
  #               stepmax = 0, hessian = T)
  return(out)
}

fitMC2 <- function(phy, x, pp = c(0.5, 0.5),
                   pi = c(0.99, 0.01), post = FALSE){

  # Test
  #sim <- TT.sim(birth = 0.2, a = 0.90, b = 0.90, tips = 100)
  #phy <- sim$tree
  #x <- sim$traits

  # Box Optimization Parameters
  min.h <- 1e-2
  max.h <- 1e2
  post <- post

  # Define Tree Features
  nb.tips <- length(phy$tip.label)
  nb.node <- phy$Nnode

  # Traits as Factors
  #x <- x[phy$tip.label]
  x <- factor(x)
  nl <- nlevels(x)
  lvls <- levels(x)
  x <- as.integer(x)

  # Define Rate Matrix
  rate <- matrix(NA, nl, nl)
  np <- nl * (nl - 1)
  rate[col(rate) != row(rate)] <- 1:np
  index.matrix <- rate
  tmp <- cbind(1:nl, 1:nl)
  index.matrix[tmp] <- NA
  rate[tmp] <- 0
  rate[rate == 0] <- np + 1

  # Define apriori rates
  pi <- pi/sum(pi)

  # Define H Matrix Output
  H <- matrix(0, nl, nl)

  # Define Likelihood Output (node -> tip)
  liks <- matrix(0, nb.tips + nb.node, nl)
  TIPS <- 1:nb.tips
  liks[cbind(TIPS, x)] <- 1

  pw <- reorder(phy, "pruningwise")

  lik <- function(pp, pi, post = FALSE){
    comp <- vector(length = nb.node + nb.tips, mode = "numeric")
    H[] <- c(pp, 0)[rate]
    diag(H) <- -rowSums(H)
    parents <- unique(pw$edge[, 1])
    root <- min(parents)
    for (i in 1:length(parents)){
      anc <- parents[i]
      ii <- which(pw$edge[, 1] == parents[i])
      desc <- pw$edge[ii, 2]
      el <- pw$edge.length[ii]
      v <- vector(length = length(desc), mode = "list")
      for (j in 1:length(v)){
        v[[j]] <- expm(H * el[j]) %*% liks[desc[j], ]}
      vv <- if(anc == root){
        Reduce("*", v)[, 1] * pi
      } else {Reduce("*", v)[, 1]}
      comp[anc] <- sum(vv)
      liks[anc, ] <- vv/comp[anc]
    }
    LogL <- -sum(log(comp[1:nb.node + nb.tips]))
    if (post != TRUE){
      return(LogL)
    } else {
      return(list(LogL = LogL, liks = liks[-TIPS,]))
    }
  }
  fit <- optim(pp, function(pp) lik(pp, pi = pi),
               method = "L-BFGS-B",
               lower = rep(0, nl), upper = rep(max.h, nl))
  lik.out <- lik(pp = fit$par, pi = pi, post = TRUE)

  if (post != TRUE){
    return(fit)
  } else {
    return(list(fit = fit, liks = lik.out))
  }
}

# Fit Function from phytools
fitMk <- function (tree, x, model = "SYM", fixedQ = NULL, ...) {
  if (hasArg(output.liks))
    output.liks <- list(...)$output.liks
  else output.liks <- FALSE
  if (hasArg(q.init))
    q.init <- list(...)$q.init
  else q.init <- length(unique(x))/sum(tree$edge.length)
  if (hasArg(opt.method))
    opt.method <- list(...)$opt.method
  else opt.method <- "nlminb"
  if (hasArg(min.q))
    min.q <- list(...)$min.q
  else min.q <- 1e-12
  N <- Ntip(tree)
  M <- tree$Nnode
  if (is.matrix(x)) {
    x <- x[tree$tip.label, ]
    m <- ncol(x)
    states <- colnames(x)
  }
  else {
    x <- to.matrix(x, sort(unique(x)))
    x <- x[tree$tip.label, ]
    m <- ncol(x)
    states <- colnames(x)
  }
  if (hasArg(pi))
    pi <- list(...)$pi
  else pi <- "equal"
  if (pi[1] == "equal")
    pi <- setNames(rep(1/m, m), states)
  else if (pi[1] == "estimated") {
    pi <- if (!is.null(fixedQ))
      statdist(fixedQ)
    else statdist(summary(fitMk(tree, x, model), quiet = TRUE)$Q)
    cat("Using pi estimated from the stationary distribution of Q assuming a flat prior.\npi =\n")
    print(round(pi, 6))
    cat("\n")
  }
  else pi <- pi/sum(pi)
  if (is.null(fixedQ)) {
    if (is.character(model)) {
      rate <- matrix(NA, m, m)
      if (model == "ER") {
        k <- rate[] <- 1
        diag(rate) <- NA
      }
      else if (model == "ARD") {
        k <- m * (m - 1)
        rate[col(rate) != row(rate)] <- 1:k
      }
      else if (model == "SYM") {
        k <- m * (m - 1)/2
        ii <- col(rate) < row(rate)
        rate[ii] <- 1:k
        rate <- t(rate)
        rate[ii] <- 1:k
      }
    }
    else {
      if (ncol(model) != nrow(model))
        stop("model is not a square matrix")
      if (ncol(model) != ncol(x))
        stop("model does not have the right number of columns")
      rate <- model
      k <- max(rate)
    }
    Q <- matrix(0, m, m)
  }
  else {
    rate <- matrix(NA, m, m)
    k <- m * (m - 1)
    rate[col(rate) != row(rate)] <- 1:k
    Q <- fixedQ
  }
  index.matrix <- rate
  tmp <- cbind(1:m, 1:m)
  rate[tmp] <- 0
  rate[rate == 0] <- k + 1
  liks <- rbind(x, matrix(0, M, m, dimnames = list(1:M + N,
                                                   states)))
  pw <- reorder(tree, "pruningwise")
  lik <- function(pp, output.liks = FALSE, pi) {
    if (any(is.nan(pp)) || any(is.infinite(pp)))
      return(1e+50)
    comp <- vector(length = N + M, mode = "numeric")
    Q[] <- c(pp, 0)[rate]
    diag(Q) <- -rowSums(Q)
    parents <- unique(pw$edge[, 1])
    root <- min(parents)
    for (i in 1:length(parents)) {
      anc <- parents[i]
      ii <- which(pw$edge[, 1] == parents[i])
      desc <- pw$edge[ii, 2]
      el <- pw$edge.length[ii]
      v <- vector(length = length(desc), mode = "list")
      for (j in 1:length(v)) v[[j]] <- matexpo(Q * el[j]) %*%
        liks[desc[j], ]
      vv <- if (anc == root)
        Reduce("*", v)[, 1] * pi
      else Reduce("*", v)[, 1]
      comp[anc] <- sum(vv)
      liks[anc, ] <- vv/comp[anc]
    }
    if (output.liks)
      return(liks[1:M + N, , drop = FALSE])
    logL <- -sum(log(comp[1:M + N]))
    return(if (is.na(logL)) Inf else logL)
  }
  if (is.null(fixedQ)) {
    if (length(q.init) != k)
      q.init <- rep(q.init[1], k)
    if (opt.method == "optim")
      fit <- optim(q.init, function(p) lik(p, pi = pi),
                   method = "L-BFGS-B", lower = rep(min.q, k))
    else fit <- nlminb(q.init, function(p) lik(p, pi = pi),
                       lower = rep(0, k), upper = rep(1e+50, k))
    obj <- list(logLik = if (opt.method == "optim") -fit$value else -fit$objective,
                rates = fit$par, index.matrix = index.matrix, states = states,
                pi = pi, method = opt.method)
    if (output.liks)
      obj$lik.anc <- lik(obj$rates, TRUE, pi = pi)
  }
  else {
    fit <- lik(Q[sapply(1:k, function(x, y) which(x == y),
                        index.matrix)], pi = pi)
    obj <- list(logLik = -fit, rates = Q[sapply(1:k, function(x,
                                                              y) which(x == y), index.matrix)], index.matrix = index.matrix,
                states = states, pi = pi)
    if (output.liks)
      obj$lik.anc <- lik(obj$rates, TRUE, pi = pi)
  }
  class(obj) <- "fitMk"
  return(obj)
}
