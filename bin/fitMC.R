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
  rate[col(rate) != row(rate)] <- np:1
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

fitMC2 <- function(phy, x, init.parms = c(0.5, 0.5),
                   prior = c(0.99, 0.01), posterior = FALSE){

  # Test
  #sim <- TT.sim(birth = 0.2, a = 0.90, b = 0.90, tips = 100)
  #phy <- sim$tree
  #x <- sim$traits

  # Box Optimization Parameters
  min.h <- 1e-5
  max.h <- 1e5
  
  # Define Parametesr
  init.parms <- init.parms
  prior <- prior
  posterior <- posterior

  # Define Tree Features
  nb.tips <- length(phy$tip.label)
  nb.node <- phy$Nnode

  # Traits as Factors
  x <- x[phy$tip.label]
  x <- factor(x)
  nl <- nlevels(x)
  lvls <- levels(x)
  x <- as.integer(x)

  # Define Rate Matrix
  rate <- matrix(NA, nl, nl)
  np <- nl * (nl - 1)
  rate[col(rate) != row(rate)] <- np:1
  index.matrix <- rate
  tmp <- cbind(1:nl, 1:nl)
  index.matrix[tmp] <- NA
  rate[tmp] <- 0
  rate[rate == 0] <- np + 1

  # Define apriori rates
  prior <- prior/sum(prior)

  # Define H Matrix Output
  H <- matrix(0, nl, nl)

  # Define Likelihood Output (node -> tip)
  liks <- matrix(0, nb.tips + nb.node, nl)
  TIPS <- 1:nb.tips
  liks[cbind(TIPS, x)] <- 1

  pw <- reorder(phy, "pruningwise")

  lik <- function(init.parms, prior, posterior = FALSE){
    if (any(is.nan(init.parms)) || any(is.infinite(init.parms)))
      return(max.h)
    comp <- vector(length = nb.node + nb.tips, mode = "numeric")
    H[] <- c(init.parms, 0)[rate]
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
        Reduce("*", v)[, 1] * prior
      } else {Reduce("*", v)[, 1]}
      comp[anc] <- sum(vv)
      liks[anc, ] <- vv/comp[anc]
    }
    LogL <- -sum(log(comp[1:nb.node + nb.tips]))
    if (posterior != TRUE){
      return(LogL)
    } else {
      return(list(LogL = LogL, liks = liks[-TIPS,]))
    }
  }
  fit <- optim(init.parms, function(init.parms) lik(init.parms, prior = prior),
               method = "L-BFGS-B",
               lower = rep(min.h, nl), upper = rep(max.h, nl))
  if (fit$convergence != 0){
    stop("Error in optim function")
  }
  lik.out <- lik(init.parms = fit$par, prior = prior, posterior = TRUE)

  if (posterior != TRUE){
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
  else min.q <- 1e-50
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
        rate[col(rate) != row(rate)] <- k:1
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
                       lower = rep(01e-50, k), upper = rep(1e+50, k))
    obj <- list(logLik = if (opt.method == "optim") -fit$value else -fit$objective,
                rates = fit$par, index.matrix = index.matrix, states = states,
                pi = pi, method = opt.method)
    if (output.liks)
      obj$lik.anc <- lik(obj$rates, TRUE, pi = pi)
  }
  else {
    fit <- lik(Q[sapply(1:k, function(x, y) which(x == y),
                        index.matrix)], pi = pi)
    obj <- list(logLik = -fit, 
                rates = Q[sapply(1:k, 
                                 function(x, y) which(x == y), 
                                 index.matrix)], index.matrix = index.matrix,
                states = states, pi = pi)
    if (output.liks)
      obj$lik.anc <- lik(obj$rates, TRUE, pi = pi)
  }
  class(obj) <- "fitMk"
  return(obj)
}

ASRTrait <- function(tree, traits, ...){
  
  # Test
  #sim <- Trait.sim(tree = tree, x = rep.comb[i,1], y = rep.comb[i,2])
  #tree <- sim$tree
  #attributes(tree)$seed <- NULL
  #traits <- sim$traits
  
  ASR <- fitMk(tree = tree, x = traits, model = "ARD", 
               output.liks = TRUE, pi = c(0.001, 0.999))
  
  tree$node.label <- as.character(1:tree$Nnode + length(tree$tip.label))
  root.dists <- as.matrix(dist.nodes(tree))[,length(tree$tip.label) + 1]

  # ID all subtrees
  subtree <- subtrees(tree, wait = FALSE)

  # Initializing Results Table
  y = rep(NA, (length(subtree)))
  cluster_size_tab <- data.frame(trait = NA, subtree = NA, node = NA,
                                 distance = NA, cluster_size = NA)

  states <- c("Off", "On")
  temp.edge <- as.data.frame(tree$edge)
  colnames(temp.edge) <- c("Anc", "Des")
  temp.edge$Len <- tree$edge.length
  temp.edge$S.Anc <- NA
  temp.edge$S.Des <- NA
  temp.edge$R.Dist <- NA
  for(j in 1:dim(ASR$lik.anc)[1]){
    temp.state <- sample(states, size = 1, prob = abs(ASR$lik.anc[j, ]))
    temp.node <- rownames(ASR$lik.anc)[j]
    temp.edge$S.Anc[which(temp.edge$Anc == temp.node)] <- temp.state
    temp.edge$S.Des[which(temp.edge$Des == temp.node)] <- temp.state
  }
  for(k in 1:dim(temp.edge)[1]){
    temp.edge$R.Dist[k] <- root.dists[temp.edge$Des[k]]
  }
  evol <- temp.edge[which(temp.edge$S.Anc == "Off" & temp.edge$S.Des == "On"), ]
  
  evols <- evol$Des    
  origins <- vector(mode = "character", length = 0)
  positives <- vector(mode = "character", length = 0)
  
  # Loop through all subtrees and determine origins
  for (j in 1:length(subtree)){
    tip_names <- subtree[[j]]$tip.label
    tree_name <- subtree[[j]]$name
    if (tree_name %in% evols){
      match_test <- match(tip_names, positives)
      if (all(is.na(match_test))){
        positives <- c(positives,tip_names)
        origins <- c(origins, tree_name)
        
        rand_tips <- sample(tip_names, size = 5, replace = T)
        cluster_dist <- distRoot(subtree[[j]], rand_tips, method = c("p"))
        cluster_size <- length(subtree[[j]]$tip.label)
        cluster_size_tab[j, ] <- c(1, j, tree_name, 
                                   mean(cluster_dist), cluster_size)
        
      } else {
        if (any(is.na(match_test))) {
          print("some NAs - something is weird")
        }
      }
    }
  }
  data.out <- cluster_size_tab[complete.cases(cluster_size_tab), ]
  return(list(ASR = ASR, Origins = data.out))
}
