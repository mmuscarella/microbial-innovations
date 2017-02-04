# Define Functions
sem <- function(x, ...){sd(x, na.rm = TRUE)/sqrt(length(na.omit(x)))}

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
  
  # Trait Evolution Plot
  #plot(y.tree, "c", FALSE, no.margin = TRUE, label.offset = 4, cex = 0.3)
  #mtext("Trait Simulation", side = 3, cex = 1, outer = F)
  #nodelabels(node.traits2$Offspring, cex = 0.25, frame = "circle", bg = n.col)
  #tiplabels(pch = 22, bg = t.col, adj = c(3, 0.5))
  
  return(min.evol)
}

TraitEvol.sim <- function(birth = 0.2, a = 0.95, b = 0.98, nsim = 100){
  replicate(nsim, TraitEvol(birth, a, b))
}

 
test.a <- seq(0.55, 0.99, 0.01)
test.b <- seq(0.55, 0.99, 0.01)

test1 <- data.frame(test.a = test.a, test.b = rep(0.90, length(test.a)),
                   sim.mean = rep(NA, length(test.a)),
                   sim.sem = rep(NA, length(test.a)))

for (i in 1:dim(test1)[1]){
  sim <- TraitEvol.sim(birth = 0.2, a = test1$test.a[i], b = test1$test.b[i], 
                       nsim = 100)
  test1$sim.mean[i] <- mean(sim)
  test1$sim.sem[i] <- sem(sim)
}

test2 <- data.frame(test.a = rep(0.90, length(test.b)), test.b = test.b,
                    sim.mean = rep(NA, length(test.b)),
                    sim.sem = rep(NA, length(test.b)))

for (i in 1:dim(test2)[1]){
  sim <- TraitEvol.sim(birth = 0.2, a = test2$test.a[i], b = test2$test.b[i], 
                       nsim = 100)
  test2$sim.mean[i] <- mean(sim)
  test2$sim.sem[i] <- sem(sim)
}



layout(matrix(1:2, 1, 2))
layout.show(2)

par(mar = c(3,3,1,1) + 0.1, oma = c(1, 1, 1, 1))

plot(test1$test.a, test1$sim.mean)
plot(test2$test.b, test2$sim.mean)

test.comb <- expand.grid(test.a, test.b)

test3 <- data.frame(test.a = test.comb[,1], test.b = test.comb[,2],
                    sim.mean = rep(NA, dim(test.comb)[1]),
                    sim.sem = rep(NA, dim(test.comb)[1]))

for (i in 1:dim(test3)[1]){
  sim <- TraitEvol.sim(birth = 0.2, a = test3$test.a[i], b = test3$test.b[i], 
                       nsim = 25)
  test3$sim.mean[i] <- mean(sim)
  test3$sim.sem[i] <- sem(sim)
}

# Trait Predictions with ASR
TraitASR <- function(tree = tree, traits = traits){
  tree <- tree
  x <- traits
  
  # Table of Observed Traits  
  # x <- Obs.Traits$Traits; names(x) <- Obs.Traits$OTU

# Use ACE function for Baysian Posterior Probabilities
ASR <- ace(x, y.tree, type = "d", model = "ARD", CI = T,
           marginal = T, ip = c(1, 0))

# Sample Posterior Probabilities for Each Node
traitNames <- c("Off","On")
prob <- 
s <- sample(traitNames, size = 1, prob = prob)



# Use Simmap Function to run Markov Chain Monte Carlo for Trait Statse on Edges
# MCMC For Time in each state on each edge
ASR.2 <- make.simmap(y.tree, x, Q = "mcmc", model = "ARD", pi = c(1,0))
edge <- ASR.2$mapped.edge / tree$edge.length
pred <- c("Off", names(unlist(lapply(ASR.2$maps, tail, n = 1L))))


# Create Trait Prediction Data Frame
pred.traits <- data.frame(Parent = traits$Parent, Offspring = traits$Offspring)
pred.traits$Trait.Pred <- NA
pred.traits$Trait.Pred[1] <- "Off"

for (i in 2:dim(pred.traits)[1]){
  po.p <- strsplit(row.names(edge), ",")[[i - 1]]
  po.o <- as.vector(unlist(pred.traits[i + 1, 1:2]), mode = "character")
  if(all.equal(po.p, po.o) == FALSE){
    stop("predicted and observed trait dataframes do not match")
  }
  pred.traits$Trait.Pred[i] <- pred[i-1]
}

pred.node <- c("Off", names(unlist(lapply(ASR.2$maps, tail, n = 1L)))[1:98])
pred.node <- gsub("On", "red", gsub("Off", "gray", pred.node))