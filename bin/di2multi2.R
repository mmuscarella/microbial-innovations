di2multi2 <- function (phy, tol = 1e-08, avg = F) {
  if (is.null(phy$edge.length)) 
    stop("the tree has no branch length")
  # ind = index of branches (edges) below threshold
  ind <- which(phy$edge.length < tol & phy$edge[, 2] > length(phy$tip.label))
  # store information about nodes that will be deleted
  node2del <- phy$edge[ind, 2]
  # n = number of branches below threshold
  n <- length(ind)
  if (!n) 
    return(phy)
  if (avg == F){
    for (i in ind){
      anc <- phy$edge[i,1]
      des <- phy$edge[i,2]
      edge.anc <- which(phy$edge[i,1] == phy$edge[,2])
      edge.des <- which(phy$edge[i,2] == phy$edge[,1])
      for (k in edge.des){
        phy$edge.length[k] <- phy$edge.length[k] + phy$edge.length[i]
        phy$edge[k,1] <- phy$edge[i,1]
      }
      phy$edge.length[i] <- NA
    }
  } else {}
  phy$edge <- phy$edge[-ind, ]
  phy$edge.length <- phy$edge.length[-ind]
  phy$Nnode <- phy$Nnode - n
  sel <- phy$edge > min(node2del)
  for (i in which(sel)) phy$edge[i] <- phy$edge[i] - sum(node2del < 
                                                           phy$edge[i])
  if (!is.null(phy$node.label)) 
    phy$node.label <- phy$node.label[-(node2del - length(phy$tip.label))]
  phy
}
  
## Example Code
test.a <- rtree(n = 50)
test.b <- rcoal(n = 50)

test <- test.b

layout(matrix(1:3), 1)
par(mar=c(0,0,0,0), oma = c(1,2,1,1))
plot(test)
mtext("original", 2)
plot(di2multi(test, 0.05))
mtext("di2multi", 2)
plot(di2multi2(test, 0.05))
mtext("di2multi+", 2)
