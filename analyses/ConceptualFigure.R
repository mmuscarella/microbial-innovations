

setwd("~/GitHub/microbial-innovations/analyses/")
source("../bin/TreeSimulationFxns.R")
source("../bin/fitMC.R")
source("../bin/ConsenTrait.R")
library("diagram")

# Plot Output
png(filename="../figures/ConceptualModel.png",
    width = 1200, height = 900, res = 96*2)

layout(matrix(c(1, 1, 2, 2, 3, 3, 3, 4), ncol = 4, byrow = T),
       heights = c(1/3, 2/3))
par(mar = (c(1, 1, 3, 1) + 0.1), oma = c(0, 0, 0, 0))


tree <- Tree.sim(tips = 20, std = 100, seed = 199)

par(mar = (c(1, 2.5, 3, 1.5) + 0.1))
plot(tree, show.tip.label = F, edge.width = 1.5,
     main = "")
mtext(side = 3, "Simulated Tree", font = 2, line = 0.5)

par(mar = (c(1, 1, 3, 1) + 0.1))
M <- matrix(c("A", "1 - A", "1 - B", "B"), 2, 2, byrow = T)
plotmat(M, pos = c(2),
        lwd = 1, box.lwd = 2,
        cex.txt = 0,
        box.size = 0.1,
        box.type = "square",
        box.prop = 0.95,
        box.col = c("gray95", "gray50"),
        arr.type = "triangle",
        arr.length = 0.2,
        arr.width = 0.2,
        arr.lwd = 1.5,
        arr.pos = 0.65, 
        endhead = T,
        self.cex = 0.4,
        self.lwd = 1.5,
        self.arrpos = c(1 * pi/2, 3 * pi/2),
        self.shifty = c(-0.09, 0.09),
        self.shiftx = c(-0.121, 0.122),
        shadow.size = 0.001,
        main = "",
        name = c("Trait\nAbsent", "Trait\nPresent"))
text(0.5, 0.35, "Gain")
text(0.5, 0.65, "Loss")
mtext(side = 3, "Trait Evolution Model", font = 2, line = 0.5)

plot(tree, show.tip.label = F, edge.width = 1.5,
     main = "")
trait <- c(rep("gray95", 7), rep("gray50", 8), 
           "gray95", "gray95", "gray50", rep("gray95", 2))
tiplabels(pch = 22, bg = trait, adj = c(3, 0.5), cex = 1.5)
#nodelabels(text = c(21:39))
nodelabels(node = c(29), adj = 0.5, pch = 21, bg = "firebrick2", cex = 2)
nodelabels(node = c(30), adj = 0.5, pch = 24, bg = "cornflower blue", cex = 2)
#nodelabels(node = c(24), adj = c(-0.5, 0.85), pch = 21, bg = "red", cex = 1.5)
#nodelabels(node = c(24), adj = c(-0.5, 0.15), pch = 24, bg = "yellow", cex = 1.5)
#text(102, 21.5, "Trait", xpd = T)
mtext(side = 3, "Trait Simulation and Inference", font = 2, line = -14, outer = T)

par(mar = (c(1, 0, 3, 1) + 0.1))
plot.new()
legend(0, 1, c("Trait Absent", "Trait Present"), pch = 22, 
       pt.bg = c("gray95", "gray50"), bty = "n")

legend(0, 0.8, c("Ancestral State\nReconstruction", "Trait\nConservation"), pch = c(21, 24),
       pt.bg = c("firebrick2", "cornflower blue"), bty = "n", 
       y.intersp = 2)

dev.off()
graphics.off()

img <- readPNG("../figures/ConceptualModel.png")
grid.raster(img)

