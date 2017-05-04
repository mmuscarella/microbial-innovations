################################################################################
#
# Trait Evolution Simulation Functions:
#       Functions used for Tree and Trait Simulations
#
# Written by: Mario Muscarella
#
# Last Update: 20170213
#
# Goals:
#       1. Write Functions to Simulate Microbial Trees and Traits
#       2. Evolve Traits
#       3. Reconstruction Evolutionary Histories
#
################################################################################

# Dependencies #################################################################
################################################################################

library("geiger")
library("MASS")
library("expm")
library("phytools")



# Standard Functions ###########################################################
###############################################################################

# Standard Error of the Mean
sem <- function(x, ...){sd(x, na.rm = TRUE)/sqrt(length(na.omit(x)))}

# Markov Chain Functions #######################################################
################################################################################

# Convert Discrete to Continuous
ab_to_xy <- function(a, b){
  x <- ((a - 1)/(a + b - 2)) * log(a + b - 1)
  y <- ((b - 1)/(a + b - 2)) * log(a + b - 1)
  #x <- (a/(a + b) * log(a + b + 1))
  #y <- (b/(a + b) * log(a + b + 1))
  return(list(x = x, y = y))}

# Convert Continuous to Discrete
xy_to_ab <- function(x, y){
  a <- (y + x*exp(-x - y))/(x + y)
  b <- (x + y*exp(-x - y))/(x + y)
  #a <- (x/(x + y) * (exp(x +y) - 1))
  #b <- (y/(x + y) * (exp(x +y) - 1))
  return(list(a = a, b = b))}

# Markov Chain Trait Evolution
MCtrait <- function(t = 5, a = "", b = "", x = "", y = ""){

  if (is.numeric(a) & is.numeric(b) == TRUE){
    Q <- matrix(c(a, 1-a, 1-b, b), 2, 2, byrow = T)

    x <- ((a - 1)/(a + b - 2)) * log(a + b - 1)
    y <- ((b - 1)/(a + b - 2)) * log(a + b - 1)

    H <- matrix(c(x, -x, -y, y), 2, 2, byrow = T)

  } else {
    H <- matrix(c(x, -x, -y, y), 2, 2, byrow = T)
  }

  H <- matrix(c(x, -x, -y, y), 2, 2, byrow = T)
  C <- eigen(H)$vectors
  D <- diag(eigen(H)$values)

  M <- expm(H * t)
  M2 <- C %*% expm(D * t) %*% ginv(C)

  return(list(M = M, M2 = M2, H = H))
}


################################################################################
# Simulation Functions #########################################################

# Tree and Trait Simulation
TT.sim <- function(birth = 0.2, a = "", b = "", x = "", y = "",
                   tips = 100, std = 100, interior = FALSE){

  # Parameters
  birth = birth; a = a; b = b; tips = tips

  # Generate Yule Tree
  y.tree <-sim.bdtree(b = birth, d = 0, stop = "taxa", n = tips, ...)
  y.tree$tip.label <- paste("OTU",
                            sprintf(paste("%0", nchar(tips)+1, "d", sep = ""),
                                    seq(1:tips)), sep = "")

  # Save a few values
  Ntips <- length(y.tree$tip.label)
  Nedges <- y.tree$Nnode + Ntips

  # Standardize Branch Lengths
  tip.dist <- round(dist.nodes(y.tree)[1:Ntips, Ntips + 1], 5)
  if (var(tip.dist) > 0.1){
    stop("Tree may not be ultrametric")
  }
  std.fac <- mean(tip.dist)
  y.tree$edge.length <- (y.tree$edge.length / std.fac) * std

  # Generate Traits Matrix
  traits <- matrix(NA, nrow = Nedges, ncol = 3)
  colnames(traits) <- c("Parent", "Offspring", "Trait")
  traits <- as.data.frame(traits)

  # Define Trait States
  traitNames <- c("Off","On")

  # Define Root Ancestor Traits
  traits[1,] <- c("None", as.character(Ntips + 1), "Off")

  # Run Trait Model Given the Tree
  for (i in 1:(Nedges - 1)){
    traits[i + 1, 1:2] <- y.tree$edge[i, ]
    t <- y.tree$edge.length[i]
    p <- y.tree$edge[i,1]
    o <- y.tree$edge[i,2]
    index.p <- which(traits$Offspring == as.character(p))
    init <- traits$Trait[index.p]
    if (init == "Off"){u <- c(1,0)} else {
      if (init == "On"){u <- c(0,1)}}
    M <- MCtrait(t = t , a = a, b = b, x = x, y = y)
    prob <- round(u %*% M$M, 6)
    s <- sample(traitNames, size = 1, prob = prob)
    traits[i + 1, 1] <- p   # Parent
    traits[i + 1, 2] <- o   # Offspring
    traits[i + 1, 3] <- s   # Trait State
  }

  if(interior == FALSE){
    tree <- y.tree
    trait <- traits$Trait[which(as.numeric(traits$Offspring) <= Ntips )]
    names(trait) <- y.tree$tip.label
    return(list(tree = tree, traits = trait))
  } else {
    if(interior == TRUE){
      tree <- y.tree
      obs.trait <- traits$Trait[which(as.numeric(traits$Offspring) <= Ntips )]
      names(obs.trait) <- y.tree$tip.label
      int.trait <- traits
      return(list(tree = tree, traits = obs.trait, interior = traits))
    } else{
      print("Please Selection T or F for Interior Traits")
    }
  }
}

# Tree and Trait Simulation
Tree.sim <- function(birth = 0.2, tips = 100, std = 100, ...){

  # Parameters
  birth = birth; tips = tips

  # Generate Yule Tree
  y.tree <-sim.bdtree(b = birth, d = 0, stop = "taxa", n = tips, ...)
  y.tree$tip.label <- paste("OTU",
                            sprintf(paste("%0", nchar(tips)+1, "d", sep = ""),
                                    seq(1:tips)), sep = "")

  # Save a few values
  Ntips <- length(y.tree$tip.label)
  Nedges <- y.tree$Nnode + Ntips

  # Standardize Branch Lengths
  tip.dist <- round(dist.nodes(y.tree)[1:Ntips, Ntips + 1], 5)
  if (var(tip.dist) > 0.1){
    stop("Tree may not be ultrametric")
  }
  std.fac <- mean(tip.dist)
  y.tree$edge.length <- (y.tree$edge.length / std.fac) * std

  return(y.tree)
}

# Tree and Trait Simulation
Trait.sim <- function(tree = tree, a , b , x , y , interior = FALSE){

  # Parameters
  #a = a; b = b; x = x; y = y

  # Define Tree
  y.tree <- tree

  # Save a few values
  Ntips <- length(y.tree$tip.label)
  Nedges <- y.tree$Nnode + Ntips

  # Generate Traits Matrix
  traits <- matrix(NA, nrow = Nedges, ncol = 3)
  colnames(traits) <- c("Parent", "Offspring", "Trait")
  traits <- as.data.frame(traits)

  # Define Trait States
  traitNames <- c("Off","On")

  # Define Root Ancestor Traits
  traits[1,] <- c("None", as.character(Ntips + 1), "Off")

  # Run Trait Model Given the Tree
  for (i in 1:(Nedges - 1)){
    traits[i + 1, 1:2] <- y.tree$edge[i, ]
    t <- y.tree$edge.length[i]
    p <- y.tree$edge[i,1]
    o <- y.tree$edge[i,2]
    index.p <- which(traits$Offspring == as.character(p))
    init <- traits$Trait[index.p]
    if (init == "Off"){u <- c(1,0)} else {
      if (init == "On"){u <- c(0,1)}}
    M <- MCtrait(t = t , x = x, y = y)
    prob <- round(u %*% M$M, 6)
    s <- sample(traitNames, size = 1, prob = prob)
    traits[i + 1, 1] <- p   # Parent
    traits[i + 1, 2] <- o   # Offspring
    traits[i + 1, 3] <- s   # Trait State
  }

  if(interior == FALSE){
    tree <- y.tree
    trait <- traits$Trait[which(as.numeric(traits$Offspring) <= Ntips )]
    names(trait) <- y.tree$tip.label
    return(list(tree = tree, traits = trait))
  } else {
    if(interior == TRUE){
      tree <- y.tree
      obs.trait <- traits$Trait[which(as.numeric(traits$Offspring) <= Ntips )]
      names(obs.trait) <- y.tree$tip.label
      int.trait <- traits
      return(list(tree = tree, traits = obs.trait, interior = traits))
    } else{
      print("Please Selection T or F for Interior Traits")
    }
  }
}

# Trait Evolution Analysis
TraitEvol <- function(tree = tree, traits.int = interior){

  # Save a few values
  Ntips <- length(tree$tip.label)
  Nedges <- tree$Nnode + Ntips

  # Define Trait States
  traitNames <- c("Off","On")

  # Extract Node and Tip States
  node.traits <- traits[which(as.numeric(traits$Offspring) > Ntips), ]
  node.traits2 <- node.traits[order(as.numeric(node.traits$Offspring)), ]
  tip.traits <- traits[which(as.numeric(traits$Offspring) <= Ntips), ]
  tip.traits2 <- tip.traits[order(as.numeric(tip.traits$Offspring)), ]

  # Define Color Vectors
  n.col <- node.traits2$Trait
  n.col <- gsub("On", "red", gsub("Off", "gray", n.col))
  t.col <- tip.traits2$Trait
  t.col <- gsub("On", "red", gsub("Off", "gray", t.col))

  # Create Observed Traits Matrix
  Obs.Traits <- data.frame(OTU = y.tree$tip.label,
                           Traits = tip.traits$Trait)

  # Add Parent Trait to Trait Matrix
  traits$P.trait <- NA
  traits$P.trait[1] <- "Off"
  for (i in 2:dim(traits)[1]){
    par.id <- which(traits$Offspring == traits$Parent[i])
    traits$P.trait[i] <- traits$Trait[par.id]
  }

  # Isolate Trait Evolution Events
  trait.evol <- traits[which(traits$P.trait == "Off" & traits$Trait == "On"), ]
  trait.evol$distance <- dist.nodes(y.tree)[trait.evol$Offspring, Ntips + 1]

  # Calculate Distance for 1st Evolution
  min.evol <- try(min(trait.evol$distance))

  return(min.evol)

  # Extract Node and Tip States
  node.traits <- traits[which(as.numeric(traits$Offspring) > Ntips), ]
  node.traits2 <- node.traits[order(as.numeric(node.traits$Offspring)), ]
  tip.traits <- traits[which(as.numeric(traits$Offspring) <= Ntips), ]
  tip.traits2 <- tip.traits[order(as.numeric(tip.traits$Offspring)), ]

  # Define Color Vectors
  n.col <- node.traits2$Trait
  n.col <- gsub("On", "red", gsub("Off", "gray", n.col))
  t.col <- tip.traits2$Trait
  t.col <- gsub("On", "red", gsub("Off", "gray", t.col))

  # Create Observed Traits Matrix
  Obs.Traits <- data.frame(OTU = y.tree$tip.label,
                           Traits = tip.traits$Trait)

  # Add Parent Trait to Trait Matrix
  traits$P.trait <- NA
  traits$P.trait[1] <- "Off"
  for (i in 2:dim(traits)[1]){
    par.id <- which(traits$Offspring == traits$Parent[i])
    traits$P.trait[i] <- traits$Trait[par.id]
  }

  # Isolate Trait Evolution Events
  trait.evol <- traits[which(traits$P.trait == "Off" & traits$Trait == "On"), ]
  trait.evol$distance <- dist.nodes(y.tree)[trait.evol$Offspring, Ntips + 1]

  # Calculate Distance for 1st Evolution
  min.evol <- try(min(trait.evol$distance))


  return(list(tree = y.tree, traits = Obs.Traits, min.evol = min.evol,
              trait.evol = trait.evol))


}

# Unlinked Tree and Trait Evolution (Yule Tree): Returns 1st Evolution
TraitEvol.old <- function(birth = 0.2, a = 0.95, b = 0.98, ...){
  # Parameters
  birth = birth; a = a; b = b

  # Generate Yule Tree
  y.tree <-sim.bdtree(b = birth, d = 0, stop = "taxa", n = 100, ...)
  y.tree$tip.label <- paste("OTU", sprintf("%05d", seq(1:100)), sep = "")

  # Save a few values
  Ntips <- length(y.tree$tip.label)
  Nedges <- y.tree$Nnode + Ntips

  # Standardize Branch Lengths
  tip.dist <- round(dist.nodes(y.tree)[1:Ntips, Ntips + 1], 5)
  if (var(tip.dist) > 0.1){
    stop("Tree may not be ultrametric")
  }
  std.fac <- mean(tip.dist)
  y.tree$edge.length <- (y.tree$edge.length / std.fac) * 100

  # Generate Traits Matrix
  traits <- matrix(NA, nrow = Nedges, ncol = 3)
  colnames(traits) <- c("Parent", "Offspring", "Trait")
  traits <- as.data.frame(traits)

  # Define Trait States
  traitNames <- c("Off","On")

  # Define Root Ancestor Traits
  traits[1,] <- c("None", as.character(Ntips + 1), "Off")

  # Run Trait Model Given the Tree
  for (i in 1:(Nedges - 1)){
    traits[i + 1, 1:2] <- y.tree$edge[i, ]
    t <- y.tree$edge.length[i]
    p <- y.tree$edge[i,1]
    o <- y.tree$edge[i,2]
    index.p <- which(traits$Offspring == as.character(p))
    init <- traits$Trait[index.p]
    if (init == "Off"){u <- c(1,0)} else {
      if (init == "On"){u <- c(0,1)}}
    M <- MCtrait(t = t , a = a, b = b)
    prob <- round(u %*% M$M, 6)
    s <- sample(traitNames, size = 1, prob = prob)
    traits[i + 1, 1] <- p   # Parent
    traits[i + 1, 2] <- o   # Offspring
    traits[i + 1, 3] <- s   # Trait State
  }

  # Above is idential to TTsim

  # Extract Node and Tip States
  node.traits <- traits[which(as.numeric(traits$Offspring) > Ntips), ]
  node.traits2 <- node.traits[order(as.numeric(node.traits$Offspring)), ]
  tip.traits <- traits[which(as.numeric(traits$Offspring) <= Ntips), ]
  tip.traits2 <- tip.traits[order(as.numeric(tip.traits$Offspring)), ]

  # Define Color Vectors
  n.col <- node.traits2$Trait
  n.col <- gsub("On", "red", gsub("Off", "gray", n.col))
  t.col <- tip.traits2$Trait
  t.col <- gsub("On", "red", gsub("Off", "gray", t.col))

  # Create Observed Traits Matrix
  Obs.Traits <- data.frame(OTU = y.tree$tip.label,
                           Traits = tip.traits$Trait)

  # Add Parent Trait to Trait Matrix
  traits$P.trait <- NA
  traits$P.trait[1] <- "Off"
  for (i in 2:dim(traits)[1]){
    par.id <- which(traits$Offspring == traits$Parent[i])
    traits$P.trait[i] <- traits$Trait[par.id]
  }

  # Isolate Trait Evolution Events
  trait.evol <- traits[which(traits$P.trait == "Off" & traits$Trait == "On"), ]
  trait.evol$distance <- dist.nodes(y.tree)[trait.evol$Offspring, Ntips + 1]

  # Calculate Distance for 1st Evolution
  min.evol <- try(min(trait.evol$distance))

  return(min.evol)
}

# Calculate First Evolution of Trait

# Unlinked Tree and Trait Evolution (Yule Tree): Returns tree and traits
TraitEvol2 <- function(birth = 0.2, a = 0.95, b = 0.98){
  # Parameters
  birth = birth; a = a; b = b

  # Generate Yule Tree
  y.tree <-sim.bdtree(b = birth, d = 0, stop = "taxa", n = 100, ...)
  y.tree$tip.label <- paste("OTU", sprintf("%05d", seq(1:100)), sep = "")

  # Save a few values
  Ntips <- length(y.tree$tip.label)
  Nedges <- y.tree$Nnode + Ntips

  # Standardize Branch Lengths
  tip.dist <- round(dist.nodes(y.tree)[1:Ntips, Ntips + 1], 5)
  if (var(tip.dist) > 0.1){
    stop("Tree may not be ultrametric")
  }
  std.fac <- mean(tip.dist)
  y.tree$edge.length <- (y.tree$edge.length / std.fac) * 100

  # Generate Traits Matrix
  traits <- matrix(NA, nrow = Nedges, ncol = 3)
  colnames(traits) <- c("Parent", "Offspring", "Trait")
  traits <- as.data.frame(traits)

  # Define Trait States
  traitNames <- c("Off","On")

  # Define Root Ancestor Traits
  traits[1,] <- c("None", "101", "Off")

  # Run Trait Model Given the Tree
  for (i in 1:(Nedges - 1)){
    traits[i + 1, 1:2] <- y.tree$edge[i, ]
    t <- y.tree$edge.length[i]
    p <- y.tree$edge[i,1]
    o <- y.tree$edge[i,2]
    index.p <- which(traits$Offspring == as.character(p))
    init <- traits$Trait[index.p]
    if (init == "Off"){u <- c(1,0)} else {
      if (init == "On"){u <- c(0,1)}}
    M <- MCtrait(t = t , a = a, b = b)
    prob <- M # round(u %*% M$M, 4)
    s <- sample(traitNames, size = 1, prob = prob)
    traits[i + 1, 1] <- p   # Parent
    traits[i + 1, 2] <- o   # Offspring
    traits[i + 1, 3] <- s   # Trait State
  }

  # Extract Node and Tip States
  node.traits <- traits[which(as.numeric(traits$Offspring) > Ntips), ]
  node.traits2 <- node.traits[order(as.numeric(node.traits$Offspring)), ]
  tip.traits <- traits[which(as.numeric(traits$Offspring) <= Ntips), ]
  tip.traits2 <- tip.traits[order(as.numeric(tip.traits$Offspring)), ]

  # Define Color Vectors
  n.col <- node.traits2$Trait
  n.col <- gsub("On", "red", gsub("Off", "gray", n.col))
  t.col <- tip.traits2$Trait
  t.col <- gsub("On", "red", gsub("Off", "gray", t.col))

  # Create Observed Traits Matrix
  Obs.Traits <- data.frame(OTU = y.tree$tip.label,
                           Traits = tip.traits$Trait)

  # Add Parent Trait to Trait Matrix
  traits$P.trait <- NA
  traits$P.trait[1] <- "Off"
  for (i in 2:dim(traits)[1]){
    par.id <- which(traits$Offspring == traits$Parent[i])
    traits$P.trait[i] <- traits$Trait[par.id]
  }

  # Isolate Trait Evolution Events
  trait.evol <- traits[which(traits$P.trait == "Off" & traits$Trait == "On"), ]
  trait.evol$distance <- dist.nodes(y.tree)[trait.evol$Offspring, Ntips + 1]

  # Calculate Distance for 1st Evolution
  min.evol <- try(min(trait.evol$distance))


  return(list(tree = y.tree, traits = Obs.Traits, min.evol = min.evol,
              trait.evol = trait.evol))
}




################################################################################
# Simulations ##################################################################

# Replicated Trait Evolution Simulations
TraitEvol.sim <- function(birth = 0.2, a = 0.95, b = 0.98, nsim = 100){
  replicate(n = nsim, expr = TraitEvol(birth, a, b))
}

# Replicated Trait Evolution Simulation with ASR
TraitEvol.sim.ASR <- function(birth = birth, a = a, b = b){
  temp <- TraitEvol2(birth, a, b)
  tree <- temp$tree
  traits <- temp$traits
  ASR <- ace(traits$Traits, tree, type = "d", model = "ARD", CI = TRUE,
         marginal = TRUE, use.expm = TRUE)
  return(ASR)
}

TraitEvol.sim.ASR2 <- function(birth = birth, a = a, b = b,
                               init.parms = c(0.9, 0.1), prior = c(0.99, 0.01)){
  temp <- TraitEvol2(birth, a, b)
  tree <- temp$tree
  traits <- temp$traits
  rownames(traits) <- traits$OTU

  ASR <- fitMC2(phy = tree, x = traits$Traits, pp = pp, pi = pi, post = TRUE)
  return(ASR)
}

# Calculate the Root State Likelihood with ACE
TraitEvol.sim2 <- function(birth = 0.2, a = 0.95, b = 0.98, nsim = 100){
  replicate(n = nsim, expr = TraitEvol.sim.ASR(birth, a, b)$lik.anc[1, ])
}


# Calculate the Root State Likelihood with fitMC2
TraitEvol.sim3 <- function(birth = 0.2, a = 0.95, b = 0.98, nsim = 100){
  replicate(n = nsim, expr = TraitEvol.sim.ASR2(birth, a, b)$liks$liks[1, ])
}

# Calculate the Posterior Likelihoods and Save Output
TraitEvolASR.Sim <- function(birth = 0.2, a = 0.95, b = 0.98, nsim = 100,
                          init.parms = c(0.5, 0.5), prior = c(0.99, 0.01)){

  SimFun <- function(birth, a, b, init.parms, prior){

    # Run Tree and Trait Simulation
    temp <- TT.sim(birth, a, b, tips = 1000)
    tree <- temp$tree
    traits <- temp$traits
    #names(traits) <- temp$traits$OTU

    # Run Ancestral State Reconstruction
    #ASR <- fitMC2(phy = tree, x = traits, init.parms = init.parms,
    #              prior = prior, posterior = TRUE)

    ASR <- fitMk(tree, traits, model = "ARD",
                 pi = prior, output.liks = TRUE)

    # Isolate Output
    pars <- -ASR$rates
    LogL <- ASR$logLik
    posterior <- ASR$lik.anc

    return(list(pars = pars, LogL = LogL, posterior = posterior))
  }
  out <- matrix(NA, nsim, 3 + 999)
  colnames(out) <- c("X", "Y", "LogL",
                     seq(1:999) + 1000)
  for(s in 1:nsim){
    temp <- try(SimFun(birth, a, b, init.parms, prior))
    if(is.list(temp)){
      out[s, 1] <- temp$pars[1]
      out[s, 2] <- temp$pars[2]
      out[s, 3] <- temp$LogL
      out[s, 4:dim(out)[2]] <- temp$posterior[1]
    }
  }
  return(out)
}

# Calculate the Trait Conservation and Save Output
TraitEvolCon.Sim <- function(birth = 0.2, a = 0.95, b = 0.98, nsim = 100,
                             level = 0.90){
  SimFun <- function(birth, a, b, level){

    # Run Tree and Trait Simulation
    temp <- TraitEvol2(birth, a, b)
    tree <- temp$tree
    attributes(tree)$seed <- NULL
    traits <- temp$traits
    rownames(traits) <- traits$OTU
    for (i in 2:dim(traits)[2]){
      traits[, i] <- as.numeric(gsub("On", 1, gsub("Off", 0, traits[, i])))
    }
    obs.first <- temp$min.evol
    obs.traits <- temp$trait.evol
    obs.numevol <- dim(temp$trait.evol)[1]

    # Run ConsenTrait
    cons <- ConsenTrait(tree = tree, traits = traits, status = FALSE)

    # Isolate Output
    cons.nodes <- cons$node
    cons.Nnodes <- dim(cons)[1]
    cons.first <- min(100 - cons$distance)

    return(c(obs.first, cons.first, obs.numevol, cons.Nnodes))
  }
  out <- matrix(NA, nsim, 4)
  colnames(out) <- c("first.obs", "first.pred", "nevol.obs", "nevol.pred")
  for (s in 1:nsim){
    out[s, ] <- try(SimFun(birth, a, b, level))
  }
  #replicate(n = nsim, expr = try(SimFun(birth, a, b, level)))
  return(out)
}
