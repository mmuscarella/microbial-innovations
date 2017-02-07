# Tree Simulation Functions

# Standard Error of the Mean
sem <- function(x, ...){sd(x, na.rm = TRUE)/sqrt(length(na.omit(x)))}

# Markov Chain Trait Evolution
MCtrait <- function(t = 5, a = 0.7, b = 0.5){
  Q <- matrix(c(a, 1-a, 1-b, b), 2, 2, byrow = T)

  x <- ((a - 1)/(a + b - 2)) * log(a + b - 1)
  y <- ((b - 1)/(a + b - 2)) * log(a + b - 1)

  H <- matrix(c(x, -x, -y, y), 2, 2, byrow = T)
  C <- eigen(H)$vectors
  D <- diag(eigen(H)$values)

  M <- expm(H * t)
  M2 <- C %*% expm(D * t) %*% ginv(C)

  return(list(M = M2, Q = Q))
}

# Unlinked Tree and Trait Evolution (Yule Tree): Returns 1st Evolution
TraitEvol <- function(birth = 0.2, a = 0.95, b = 0.98){
  # Parameters
  birth = birth; a = a; b = b

  # Generate Yule Tree
  y.tree <-sim.bdtree(b = birth, d = 0, stop = "taxa", n = 100)
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
    prob <- round(u %*% M$M, 4)
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
  min.evol <- min(trait.evol$distance)

  return(min.evol)
}

# Replicated Trait Evolution Simulations
TraitEvol.sim <- function(birth = 0.2, a = 0.95, b = 0.98, nsim = 100){
  replicate(nsim, TraitEvol(birth, a, b))
}


# Unlinked Tree and Trait Evolution (Yule Tree): Returns tree and traits
TraitEvol2 <- function(birth = 0.2, a = 0.95, b = 0.98){
  # Parameters
  birth = birth; a = a; b = b

  # Generate Yule Tree
  y.tree <-sim.bdtree(b = birth, d = 0, stop = "taxa", n = 100)
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
    prob <- round(u %*% M$M, 4)
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
  min.evol <- min(trait.evol$distance)


  return(list(tree = y.tree, traits = Obs.Traits, min.evol = min.evol))
}

TraitEvol.sim.ASR <- function(birth = birth, a = a, b = b){
  temp <- TraitEvol2(birth, a, b)
  tree <- temp$tree
  traits <- temp$traits
  ASR <- ace(traits$Traits, tree, type = "d", model = "ARD", CI = TRUE,
         marginal = TRUE, use.expm = TRUE)
  return(ASR)
}

# Calculate the Root State Likelihood
TraitEvol.sim2 <- function(birth = 0.2, a = 0.95, b = 0.98, nsim = 100){
  replicate(nsim, TraitEvol.sim.ASR(birth, a, b)$lik.anc[1, ])
}



# Tree and Trait Simulation
TT.sim <- function(birth = 0.2, a = 0.90, b = 0.90, tips = 100){
  # Parameters
  birth = birth; a = a; b = b; tips = tips
  
  # Generate Yule Tree
  y.tree <-sim.bdtree(b = birth, d = 0, stop = "taxa", n = tips)
  y.tree$tip.label <- paste("OTU", sprintf("%05d", seq(1:tips)), sep = "")
  
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
    prob <- round(u %*% M$M, 4)
    s <- sample(traitNames, size = 1, prob = prob)
    traits[i + 1, 1] <- p   # Parent
    traits[i + 1, 2] <- o   # Offspring
    traits[i + 1, 3] <- s   # Trait State
  }
  tree <- y.tree
  trait <- traits$Trait[which(as.numeric(traits$Offspring) <= tips )]
  names(trait) <- y.tree$tip.label
  
  return(list(tree = tree, traits = trait))

}


