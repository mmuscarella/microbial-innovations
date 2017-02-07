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
  liks <- matrix(0, nb.tip + nb.node, nl)
  TIPS <- 1:nb.tip
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
    comp <- numeric(nb.tip + nb.node)
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

fitMC2 <- function(phy, x, pp = 0.5, pi = c(1, 0)){
  
  # Test 
  #sim <- TT.sim(birth = 0.2, a = 0.90, b = 0.90, tips = 100)
  #phy <- sim$tree
  #x <- sim$traits
  
  # Box Optimization Parameters
  min.h <- 1e-15
  max.h <- 10e15
  
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
  liks <- matrix(0, nb.tip + nb.node, nl)
  TIPS <- 1:nb.tips
  liks[cbind(TIPS, x)] <- 1

  pw <- reorder(phy, "pruningwise")
  
  lik <- function(pp, pi){
    comp <- vector(length = nb.node + nb.tip, mode = "numeric")
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
    return(LogL)
  }
    
  fit <- optim(rep(pp, length.out = nl), function(pp) lik(pp, pi = pi), 
               method = "L-BFGS-B", 
               lower = rep(min.h, nl), upper = rep(max.h, nl))
  return(fit)
}

a <- 0.11
b <- 0.90
sim <- TT.sim(birth = 0.2, a = a, b = b, tips = 100)
Hx <- ((a - 1)/(a + b - 2)) * log(a + b - 1)
Hy <- ((b - 1)/(a + b - 2)) * log(a + b - 1)
phy <- sim$tree
x <- sim$traits

fitMC(phy = phy, x = x, ip = 0.5)
fitMC2(phy = phy, x = x, pp = 0.5, pi = c(1, 0))
  
Hx; Hy

fitMC <- function(x, phy, ip = 0.5){
  
  # Test 
  sim <- TT.sim(birth = 0.2, a = 0.90, b = 0.90, tips = 100)
  phy <- sim$tree
  x <- sim$traits
  ip = 0.5
  
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
  liks <- matrix(0, nb.tip + nb.node, nl)
  TIPS <- 1:nb.tip
  liks[cbind(TIPS, x)] <- 1
  
  phy <- reorder(phy, "postorder")
  H <- matrix(0, nl, nl)
  e1 <- phy$edge[, 1]
  e2 <- phy$edge[, 2]
  EL <- phy$edge.length
  
  # Define Function
  E <- expm::expm
  lik <- function(p){
    if (any(is.nan(p)) || any(is.infinite(p)))
      return(1e+50)
    comp <- numeric(nb.tip + nb.node)
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
    l <- -sum(log(comp[-TIPS]))
    
    return(l)
  }

