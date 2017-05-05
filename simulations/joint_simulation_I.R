#rm(list=ls())
#setwd("~/GitHub/microbial-innovations/simulations")
setwd("~/microbial-innovations/simulations")
source("../bin/TreeSimulationFxns.R")
source("../bin/fitMC.R")
source("../bin/ConsenTrait.R")
library("methods")
library("adephylo")

# Add Basic Functions
se <- function(x, ...){sd(x, na.rm = TRUE)/sqrt(length(na.omit(x)))}
ci <- function(x, ...){1.96 * sd(x,na.rm = TRUE)}


# Simulation 1: Accuracy of ASR Predictions
test.x <- c(-0.2, -0.02, -0.002, -0.0002, -0.00002, -0.000002)
test.y <- c(-0.2, -0.02, -0.002, -0.0002, -0.00002, -0.000002)
test.comb <- data.frame(expand.grid(X = test.x, Y = test.y))
rep.comb <- data.frame(test.comb[rep(seq_len(nrow(test.comb)), each=10),])
rownames(rep.comb) <- NULL

mc.testA <- matrix(NA, ncol = 6, nrow = dim(rep.comb)[1])
colnames(mc.testA) <- c("X", "Y", "X.P", "Y.P", "L.On", "L.Off")

tree <- Tree.sim(tips = 2000, std = 4000)

for (i in 1:dim(mc.testA)[1]){
  print(paste("Simulation", i, "out of", dim(mc.testA)[1], sep = " "))
  sim <- Trait.sim(tree = tree, x = rep.comb[i,1], y = rep.comb[i,2])
  tree <- sim$tree
  traits <- sim$traits
  mc.testA[i, 1] <- rep.comb[i,1]
  mc.testA[i, 2] <- rep.comb[i,2]
  mc.testA[i, 5] <- sum(traits == "On")
  mc.testA[i, 6] <- sum(traits == "Off")
  
  if(sum(traits == "On") > 0 & sum(traits == "On") < length(traits)){
  # Run Ancestral State Reconstruction
  ASR <- fitMk(tree, traits, model = "ARD", 
               output.liks = TRUE, pi = c(0.001, 0.999))
  
  # Isolate Output
  rates <- -ASR$rates
  mc.testA[i, 3:4] <- c(rates[1], rates[2])
  } else {
    mc.testA[i, 3:4] <- c(NA, NA)
  }
}



mc.testA[which(mc.testA[, 3] < -0.000001), ]
mc.testA2 <- mc.testA[which(mc.testA[, 5] > 2 & mc.testA[, 6] > 2), ]

# Plot Output
library("grid")
library("png")
png(filename="../simulations/output/Estimates.png",
    width = 1800, height = 900, res = 96*2)

layout(matrix(1:2, ncol = 2, byrow = T))
par(mar = c(4, 4, 3, 1) + 0.5, oma = c(0, 0, 0, 0))

plot(y = log10(-mc.testA[, 3]), x = jitter(log10(-mc.testA[, 1])),
     xlab = "", ylab = "",
     xlim = c(-6.5, 0), ylim = c(-6.5, 0), axes = F,
     col=rgb(0, 0, 1, 0.5), pch=16)

mtext(side = 1, expression(paste("True Rate (log"[10], ")")), 
      line = 2.5, cex = 1.2)
mtext(side = 2, expression(paste("Predicted Rate (log"[10], ")")), 
      line = 2.5, cex = 1.2)
mtext(side = 3, text = "X Rate", line = 1, cex = 1.25)

abline(0, 1, lty = 2)

axis(side = 1, labels = T, at = c(-6, -4, -2, 0), las = 1, 
     lwd.ticks = 2, tck = -0.02, lwd = 1.5)
axis(side = 2, labels = T, at = c(-6, -4, -2, 0), las = 1, 
     lwd.ticks = 2, tck = -0.02, lwd = 1.5)
axis(side = 3, labels = F, at = c(-6, -4, -2, 0), las = 1, 
     lwd.ticks = 2, tck = -0.02, lwd = 1.5)
axis(side = 4, labels = F, at = c(-6, -4, -2, 0), las = 1, 
     lwd.ticks = 2, tck = -0.02, lwd = 1.5)
axis(side = 1, labels = F, at = c(-6, -4, -2, 0), las = 1, 
     lwd.ticks = 2, tck = 0.02, lwd = 1.5)
axis(side = 2, labels = F, at = c(-6, -4, -2, 0), las = 1, 
     lwd.ticks = 2, tck = 0.02, lwd = 1.5)
axis(side = 3, labels = F, at = c(-6, -4, -2, 0), las = 1, 
     lwd.ticks = 2, tck = 0.02, lwd = 1.5)
axis(side = 4, labels = F, at = c(-6, -4, -2, 0), las = 1, 
     lwd.ticks = 2, tck = 0.02, lwd = 1.5)
box(lwd = 1.5)

plot(y = log10(-mc.testA[, 4]), x = jitter(log10(-mc.testA[, 2])),
     xlab = "", ylab = "",
     xlim = c(-6.5, 0), ylim = c(-6.5, 0), axes = F,
     col=rgb(0, 0, 1, 0.5), pch=16)

mtext(side = 1, expression(paste("True Rate (log"[10], ")")), 
      line = 2.5, cex = 1.2)
mtext(side = 2, expression(paste("Predicted Rate (log"[10], ")")), 
      line = 2.5, cex = 1.2)
mtext(side = 3, text = "Y Rate", line = 1, cex = 1.25)

abline(0, 1, lty = 2)

axis(side = 1, labels = T, at = c(-6, -4, -2, 0), las = 1, 
     lwd.ticks = 2, tck = -0.02, lwd = 1.5)
axis(side = 2, labels = T, at = c(-6, -4, -2, 0), las = 1, 
     lwd.ticks = 2, tck = -0.02, lwd = 1.5)
axis(side = 3, labels = F, at = c(-6, -4, -2, 0), las = 1, 
     lwd.ticks = 2, tck = -0.02, lwd = 1.5)
axis(side = 4, labels = F, at = c(-6, -4, -2, 0), las = 1, 
     lwd.ticks = 2, tck = -0.02, lwd = 1.5)
axis(side = 1, labels = F, at = c(-6, -4, -2, 0), las = 1, 
     lwd.ticks = 2, tck = 0.02, lwd = 1.5)
axis(side = 2, labels = F, at = c(-6, -4, -2, 0), las = 1, 
     lwd.ticks = 2, tck = 0.02, lwd = 1.5)
axis(side = 3, labels = F, at = c(-6, -4, -2, 0), las = 1, 
     lwd.ticks = 2, tck = 0.02, lwd = 1.5)
axis(side = 4, labels = F, at = c(-6, -4, -2, 0), las = 1, 
     lwd.ticks = 2, tck = 0.02, lwd = 1.5)
box(lwd = 1.5)


dev.off()
graphics.off()

img <- readPNG("../simulations/output/Estimates.png")
grid.raster(img)



# Simulation 2: How does 1st Evolution Change with X
test.x <- c(-0.2, -0.02, -0.002, -0.0002, -0.00002, -0.000002)
test.y <- c(-0.2, -0.02, -0.002, -0.0002, -0.00002, -0.000002)
test.comb <- data.frame(expand.grid(X = test.x, Y = test.y))
rep.comb <- data.frame(test.comb[rep(seq_len(nrow(test.comb)), each=10),])
rownames(rep.comb) <- NULL

mc.testB <- data.frame(matrix(NA, ncol = 12, nrow = dim(rep.comb)[1]))
colnames(mc.testB) <- c("X", "Y", "X.P", "Y.P", "L.On", "L.Off", 
                       "ASR.Origins", "ASR.Mean", "ASR.SE",
                       "Con.Origins", "Con.Mean", "Con.SE")

tree <- Tree.sim(tips = 2000, std = 4000)

for (i in 1:dim(mc.testB)[1]){
  print(paste("Simulation", i, "out of", dim(mc.testB)[1], sep = " "))
  sim <- Trait.sim(tree = tree, x = rep.comb[i,1], y = rep.comb[i,2])
  tree <- sim$tree
  attributes(tree)$seed <- NULL
  traits <- sim$traits
  mc.testB[i, 1] <- rep.comb[i,1]
  mc.testB[i, 2] <- rep.comb[i,2]
  mc.testB[i, 5] <- sum(traits == "On")
  mc.testB[i, 6] <- sum(traits == "Off")
  
  if(sum(traits == "On") > 2 & sum(traits == "Off") > 2){
    # Run Ancestral State Reconstruction
    ASR <- ASRTrait(tree, traits)
    ASR.r <- ASR$ASR
    ASR.o <- ASR$Origins
    
    # Run ConsenTrait
    trait.tab <- data.frame(names(traits), (traits == "On"))
    CON <- ConsenTrait(tree = tree, traits = trait.tab)
    
    # Isolate Output
    rates <- -ASR.r$rates
    mc.testB[i, 3:4] <- c(rates[1], rates[2])
    mc.testB[i, 7] <- try(dim(ASR.o)[1])
    mc.testB[i, 8] <- try(mean(ASR.o$distance))
    mc.testB[i, 9] <- try(se(ASR.o$distance))
    mc.testB[i, 10] <- try(dim(CON)[1])
    mc.testB[i, 11] <- try(mean(CON$distance))
    mc.testB[i, 12] <- try(se(CON$distance))
    } else {
    mc.testB[i, c(3:4, 7:12)] <- NA
  }
}


save.image(file = "joint_simulation_I.RData")



