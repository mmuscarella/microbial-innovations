source("../bin/TreeSimulationFxns.R")

# Simulation 1: How does 1st Evolution Change with A
test.a <- seq(0.5, 0.99, 0.02)
test.b <- seq(0.5, 0.99, 0.02)

test1 <- data.frame(test.a = test.a, test.b = rep(0.90, length(test.a)),
                   sim.mean = rep(NA, length(test.a)),
                   sim.sem = rep(NA, length(test.a)))

for (i in 1:dim(test1)[1]){
  sim <- TraitEvol.sim(birth = 0.2, a = test1$test.a[i], b = test1$test.b[i], 
                       nsim = 10)
  test1$sim.mean[i] <- mean(sim)
  test1$sim.sem[i] <- sem(sim)
}

# Simulation 2: How does 1st Evolution Change with B
test2 <- data.frame(test.a = rep(0.90, length(test.b)), test.b = test.b,
                    sim.mean = rep(NA, length(test.b)),
                    sim.sem = rep(NA, length(test.b)))

for (i in 1:dim(test2)[1]){
  sim <- TraitEvol.sim(birth = 0.2, a = test2$test.a[i], b = test2$test.b[i], 
                       nsim = 10)
  test2$sim.mean[i] <- mean(sim)
  test2$sim.sem[i] <- sem(sim)
}

# Plot with Simulation 1 and 2
layout(matrix(1:2, 1, 2))
layout.show(2)

par(mar = c(4,1,1,1) + 0.5, oma = c(1, 4, 0.5, 0.5))

plot(test1$test.a, test1$sim.mean, las = 1, pch = 16,
     xlim = c(0.4, 1.0), ylim = c(5, 54),
     xlab = "Alpha Parameter", ylab = "")
arrows(test1$test.a, y0 = test1$sim.mean - test1$sim.sem,
       y1 = test1$sim.mean + test1$sim.sem, length = 0)
legend("topleft", legend = expression(paste(beta, " = 0.90")), 
       bty = "n", xjust = 0, lty = 0, adj = c(0.3, 0.5))
mtext(side = 2, "Mean 1st Occurrence (time)", cex = 1, line = 3)

plot(test2$test.b, test2$sim.mean, las = 1, pch = 16,
     xlim = c(0.4, 1.0), ylim = c(5, 54),
     xlab = "Beta Parameter", ylab = "")
arrows(test2$test.b, y0 = test2$sim.mean - test2$sim.sem,
       y1 = test2$sim.mean + test2$sim.sem, length = 0)
legend("topright", legend = expression(paste(alpha, " = 0.90")), 
       bty = "n", xjust = 0, lty = 0)

# Test 3: 1st Evolution across ranges of A and B
test.comb <- expand.grid(test.a, test.b)[-1,]

test3 <- data.frame(test.a = test.comb[,1], test.b = test.comb[,2],
                    sim.mean = rep(NA, dim(test.comb)[1]),
                    sim.sem = rep(NA, dim(test.comb)[1]))

for (i in 1:dim(test3)[1]){
  sim <- TraitEvol.sim(birth = 0.2, a = test3$test.a[i], b = test3$test.b[i], 
                       nsim = 10)
  test3$sim.mean[i] <- mean(sim)
  test3$sim.sem[i] <- sem(sim)
}

# Plot Test #3
layout(1)
par(mar = c(1, 1, 1, 1), oma = c(0, 0, 0, 0))

l <- length(unique(test3[,1]))
w <- length(unique(test3[,2]))
z <- matrix(c(NA, test3[, 3]), nrow = length(test.a), ncol = length(test.b))
nrz <- nrow(z)
ncz <- ncol(z)
# Create a function interpolating colors in the range of specified colors
jet.colors <- colorRampPalette( c("blue", "green") )
# Generate the desired number of colors from this palette
nbcol <- 100
color <- jet.colors(nbcol)
# Compute the z-value at the facet centres
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
# Recode facet z-values into color indices
facetcol <- cut(zfacet, nbcol)
par(mar = c(2, 2, 2, 2) + 0.5)
persp(z = z, phi = 15, theta = 300,
      zlim = c(0, 75), ticktype = "simple", nticks = 4,
      xlab = "Alpha", ylab = "Beta", zlab = "Mean 1st Occurence",
      col = color[facetcol], scale = TRUE)





# How does the likelihood of the root state change with A and B
test.a <- seq(0.54, 0.95, 0.05)
test.b <- seq(0.54, 0.95, 0.05)
test.comb <- expand.grid(test.a, test.b)

sim.lik <- data.frame(test.a = test.comb[,1], test.b = test.comb[,2],
                               sim.mean = rep(NA, dim(test.comb)[1]),
                               sim.sem = rep(NA, dim(test.comb)[1]))

for (i in 1:dim(test.comb)[1]){
  sim <- TraitEvol.sim2(birth = 0.2, a = test.comb[i, 1], b = test.comb[i, 2], 
                       nsim = 10)
  sim.lik$sim.mean[i] <- mean(sim[1, ])
  #sim.lik$sim.sem[i] <- sem(sim[1, ])
}
x <- test.a
y <- test.b
z <- matrix(sim.lik[, 3], nrow = length(test.a), ncol = length(test.b))
nrz <- nrow(z)
ncz <- ncol(z)
# Create a function interpolating colors in the range of specified colors
jet.colors <- colorRampPalette( c("blue", "green") )
# Generate the desired number of colors from this palette
nbcol <- 100
color <- jet.colors(nbcol)
# Compute the z-value at the facet centres
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
# Recode facet z-values into color indices
facetcol <- cut(zfacet, nbcol)
par(mar = c(2, 2, 2, 2) + 0.5)
persp(z = z, phi = 15, theta = 30,
      zlim = c(0.25, 0.8), ticktype = "simple", nticks = 4,
      xlab = "Alpha", ylab = "Beta", zlab = "Likelihood of State",
      col = color[facetcol], scale = FALSE)


# How does the likelihood of the root state change with A and B (new function)
test.a <- seq(0.54, 0.95, 0.05)
test.b <- seq(0.54, 0.95, 0.05)
test.comb <- expand.grid(test.a, test.b)

sim.lik <- data.frame(test.a = test.comb[,1], test.b = test.comb[,2],
                      sim.mean = rep(NA, dim(test.comb)[1]),
                      sim.sem = rep(NA, dim(test.comb)[1]))

for (i in 1:dim(test.comb)[1]){
  sim <- TraitEvol.sim3(birth = 0.2, a = test.comb[i, 1], b = test.comb[i, 2], 
                        nsim = 10)
  sim.lik$sim.mean[i] <- mean(sim[1, ])
  #sim.lik$sim.sem[i] <- sem(sim[1, ])
}
x <- test.a
y <- test.b
z <- matrix(sim.lik[, 3], nrow = length(test.a), ncol = length(test.b))
nrz <- nrow(z)
ncz <- ncol(z)
# Create a function interpolating colors in the range of specified colors
jet.colors <- colorRampPalette( c("blue", "green") )
# Generate the desired number of colors from this palette
nbcol <- 100
color <- jet.colors(nbcol)
# Compute the z-value at the facet centres
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
# Recode facet z-values into color indices
facetcol <- cut(zfacet, nbcol)
par(mar = c(2, 2, 2, 2) + 0.5)
persp(z = z, phi = 15, theta = 30,
      zlim = c(0.9, 1), ticktype = "simple", nticks = 4,
      xlab = "Alpha", ylab = "Beta", zlab = "Likelihood of State",
      col = color[facetcol], scale = FALSE)

# How Does the Likelihood change with A and B (new function)
test.a <- seq(0.55, 0.95, 0.05)
test.b <- seq(0.55, 0.95, 0.05)
test.comb <- expand.grid(test.a, test.b)

test.combB <- test.comb
test.combB$x <- NA; test.combB$y <- NA
for (i in 1:dim(test.combB)[1]){
  test.combB[i, 3] <- ab_to_xy(test.combB[i, 1], test.combB[i, 2])[1]
  test.combB[i, 4] <- ab_to_xy(test.combB[i, 1], test.combB[i, 2])[2]
}

sim.lik <- data.frame(test.a = test.comb[,1], test.b = test.comb[,2],
                      ll.mean = rep(NA, dim(test.comb)[1]),
                      ll.sem = rep(NA, dim(test.comb)[1]),
                      x.mean = rep(NA, dim(test.comb)[1]),
                      y.mean = rep(NA, dim(test.comb)[1]))

for (i in 1:dim(test.comb)[1]){
  sim <- TraitEvol.sim4(birth = 0.2, a = test.comb[i, 1], b = test.comb[i, 2], 
                        nsim = 10)
  sim.lik$ll.mean[i] <- mean(sapply(sim[2,], '['))
  sim.lik$ll.sem[i] <- sem(sapply(sim[2,], '['))
  sim.lik$x.mean[i] <- mean(sapply(sim[1,], '[')[1,])
  sim.lik$y.mean[i] <- mean(sapply(sim[1,], '[')[2,])
}

x <- test.a
y <- test.b
z <- matrix(sim.lik$ll.mean, nrow = length(test.a), ncol = length(test.b))


nrz <- nrow(z)
ncz <- ncol(z)
# Create a function interpolating colors in the range of specified colors
jet.colors <- colorRampPalette( c("blue", "green") )
# Generate the desired number of colors from this palette
nbcol <- 100
color <- jet.colors(nbcol)
# Compute the z-value at the facet centres
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
# Recode facet z-values into color indices
facetcol <- cut(zfacet, nbcol)
par(mar = c(2, 2, 2, 2) + 0.5)
persp(z = z, phi = 15, theta = 30,
      zlim = c(0.9, 1), ticktype = "simple", nticks = 4,
      xlab = "Alpha", ylab = "Beta", zlab = "Likelihood of State",
      col = color[facetcol], scale = FALSE)











# Trait Predictions with ASR
TraitASR <- function(tree = tree, traits = traits){
  tree <- tree
  x <- traits
  
  # Table of Observed Traits  
  # x <- Obs.Traits$Traits; names(x) <- Obs.Traits$OTU

# Use ACE function for Baysian Posterior Probabilities
ASR <- ace(x, tree, type = "d", model = "ARD", CI = T,
           marginal = T, ip = c(1, 0))

ASR.M <- matrix(c(-ASR$rates[2], ASR$rates[2], 
                     ASR$rates[1], -ASR$rates[1]), nrow = 2, ncol = 2, byrow = T)
post <- ASR$lik.anc
pred.state <- as.data.frame(matrix(NA, dim(post)[1]))
                            
p <- expm::expm(Q)                        


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