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

di2multi3 <- function (phy1 = "", phy2 = "", tol = 1e-08, avg = F) {
  if (is.null(phy1$edge.length)) 
    stop("the primary tree has no branch length")
  if (is.null(phy2$edge.length))
    stop("the secondary tree has no branch lengh")
  if (all.equal(phy1$node.label, phy2$node.label) == FALSE)
    stop("the trees must have the same nodes")
  # identify standing polytomies
  polys.size <- table(phy2$edge[,1])
  polys.old <- names(polys.size[which(polys.size > 2)])
  # ind = index of branches (edges) in primary tree below threshold
  ind <- which(phy1$edge.length < tol & phy1$edge[, 2] > length(phy1$tip.label))
  # store information about nodes that will be deleted
  node2del <- phy1$edge[ind, 2]
  # n = number of branches below threshold
  n <- length(ind)
  if (!n) 
    return(phy2)
  if (avg == F){
    for (i in ind){
      anc <- phy2$edge[i,1]
      des <- phy2$edge[i,2]
      edge.anc <- which(phy2$edge[i,1] == phy2$edge[,2])
      edge.des <- which(phy2$edge[i,2] == phy2$edge[,1])
      for (k in edge.des){
        phy2$edge.length[k] <- phy2$edge.length[k] + phy2$edge.length[i]
        phy2$edge[k,1] <- phy2$edge[i,1]
      }
      phy2$edge.length[i] <- NA
    }
  } else {
    stop("no method to average yet")
  }
  
  # Clean Up Output
  phy2$edge <- phy2$edge[-ind, ]
  phy2$edge.length <- phy2$edge.length[-ind]
  phy2$Nnode <- phy2$Nnode - n
  sel <- phy2$edge > min(node2del)
  for (i in which(sel)) phy2$edge[i] <- phy2$edge[i] - sum(node2del < 
                                                           phy2$edge[i])
  if (!is.null(phy2$node.label)) 
    phy2$node.label <- phy2$node.label[-(node2del - length(phy2$tip.label))]
  
  # identify new polytomies
  polys.size <- table(phy2$edge[,1])
  polys.temp <- names(polys.size[which(polys.size > 2)])
  polys.new <- setdiff(polys.temp, polys.old)
  
  phy2$multi.old = as.numeric(polys.old)
  phy2$multi.new = as.numeric(polys.new)
  return(phy2)
}
  
## Example Code
test.di2multi2 <- function(n = 50){
  test.a <- rtree(n = n)
  test.b <- rcoal(n = n)
  
  test <- test.b
  
  layout(matrix(1:3), 1)
  par(mar=c(0,0,0,0), oma = c(1,2,1,1))
  plot(test)
  mtext("original", 2)
  plot(di2multi(test, 0.05))
  mtext("di2multi", 2)
  plot(di2multi2(test, 0.05))
  mtext("di2multi+", 2)
}
