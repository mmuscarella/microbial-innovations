# Tree Simulations

source("../bin/TreeSimulationFnxs.R")
 
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










TraitEvol.sim2 <- function(birth = 0.2, a = 0.95, b = 0.98, nsim = 100){
  replicate(nsim, TraitEvol(birth, a, b))
}


test.a <- seq(0.55, 0.99, 0.05)
test.b <- seq(0.55, 0.99, 0.05)
test.comb <- expand.grid(test.a, test.b)

test3 <- data.frame(test.a = test.comb[,1], test.b = test.comb[,2],
                    sim.mean = rep(NA, dim(test.comb)[1]),
                    sim.sem = rep(NA, dim(test.comb)[1]))

for (i in 1:dim(test3)[1]){
  sim <- TraitEvol.sim2(birth = 0.2, a = test3$test.a[i], b = test3$test.b[i], 
                       nsim = 1)
  test3$sim.mean[i] <- mean(sim)
  test3$sim.sem[i] <- sem(sim)
}
