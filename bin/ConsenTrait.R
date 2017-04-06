ConsenTrait <- function(tree = "", traits = "", cutoff = 0.9,
                        status = TRUE){

  require("adephylo")||install.packages("adephylo");require("adephylo")
  require("phytools")||install.packages("phytools");require("phytools")
  require("ape")||install.packages("ape");require("ape")

  # Import Tree and Root  if Needed
  if (is.rooted(tree) == FALSE){
    root_tree <- midpoint.root(tree)
  } else {
    root_tree <- tree
  }

  # Import Traits into Function
  table <- traits

  # Drop tips not found in trait table
  z <- subset(tree$tip.label, !(tree$tip.label %in% table[,1]) )
  if (length(z) > 0){
    drop.tip(tree,z)
  }

  # Replace any negative branch lengths
  root_tree$edge.length[root_tree$edge.length <= 0] = 0.00001
  
  root_tree$node.label.old <- root_tree$node.label
  root_tree$node.label <- as.character(1:root_tree$Nnode + 
                                         length(root_tree$tip.label))

  # ID all subtrees
  subtree <- subtrees(root_tree, wait = FALSE)

  # Initializing Results Table
  y = rep(NA, (length(subtree) * (dim(table)[2] - 1)))
  cluster_size_tab <- data.frame(trait = NA, subtree = NA, node = NA,
                                 distance = NA, cluster_size = NA)

  # Loop Through Traits
  for (i in 2:ncol(table)){

    # Status Indicator
    if (status == TRUE){
      print(paste("Analyzing Trait", i - 1, "of",
                  ncol(table)[[1]] - 1, "...", sep = " "), quote = F)
    }

    # Make Temp Table
    table_tmp <- data.frame(ID = table[, 1], Trait = table[,i])

    # Remove All Entries Not in Tree and Sort by ID
    table2 <- table_tmp[which(table_tmp$ID %in% root_tree$tip.label), ]
    table2 <- table2[sort(table2$ID), ]

    # Initialize Temp Result Vectors
    positives <- vector(mode = "character", length = 0)
    cluster_size <- numeric(length=0)
    cluster_dist <- numeric(length = 0)
    node_positive <- vector(mode = "character", length = 0)

    # Loop through all subtrees and determine if any subtrees have >90% positives
    for (j in 1:length(subtree)){
      tip_names <- subtree[[j]]$tip.label
      if (mean(table2$Trait[which(table2$ID %in% tip_names)]) > cutoff){
        match_test <- match(tip_names, positives)
        if (all(is.na(match_test))){
          positives <- c(positives,tip_names)
          node_positive <- subtree[[j]]$node.label[1]
          
          rand_tips <- sample(tip_names, size = 5, replace = T)
          cluster_dist <- distRoot(subtree[[j]], rand_tips, method = c("p"))
          cluster_size <- length(subtree[[j]]$tip.label)
          cluster_size_tab[j + length(subtree) * (i - 2), ] <- c(i - 1, j, 
                               node_positive, mean(cluster_dist), cluster_size)

        } else {
          if (any(is.na(match_test))) {
            print("some NAs - something is weird")
          }
        }
      }
    }
  }
  data.out <- cluster_size_tab[complete.cases(cluster_size_tab), ]
  return(data.out)
}
