# https://doi.org/10.1214/22-AOAS1643
# https://doi.org/10.1214/22-AOAS1643SUPPB

## The auxiliary functions for building spanning tree

################################################################################
################################################################################
###################### For spanning-tree operations

#### Generate clusters for the simulation data

GenerateClust <- function(coord, clustsize) {
  ns <- nrow(coord)
  graph0 <- ConstructGraph0(coord, method = "knn", 10)
  E(graph0)$weights <- runif(length(E(graph0)))
  mstgraph <- mst(graph0)

  ### For balanced clusters
  bsize <- floor(ns / clustsize / 2)
  membership <- 1:ns
  while (min(table(membership)) < bsize) {
    delid <- sample(1:(vcount(mstgraph) - 1), clustsize - 1)
    mstgraph2 <- delete.edges(mstgraph, delid)
    membership <- components(mstgraph2)$membership
  }

  return(list(graph0 = graph0, truemst = mstgraph, membership = components(mstgraph2)$membership))
}

### Unbalanced clusters

GenerateClust_unbalanced <- function(coord, clustsize) {
  ns <- nrow(coord)
  graph0 <- ConstructGraph0(coord, method = "knn", 10)
  E(graph0)$weights <- runif(length(E(graph0)))
  mstgraph <- mst(graph0)

  ### For balanced clusters
  ### bsize=floor(ns/clustsize/2)
  bsize <- 6
  membership <- 1:ns
  while (min(table(membership)) > bsize | min(table(membership)) < 2) {
    delid <- sample(1:(vcount(mstgraph) - 1), clustsize - 1)
    mstgraph2 <- delete.edges(mstgraph, delid)
    membership <- components(mstgraph2)$membership
  }

  return(list(graph0 = graph0, truemst = mstgraph, membership = components(mstgraph2)$membership))
}

# function to get whether an edge is within a cluster or bewteen two clusters
getEdgeStatus <- function(membership, graph) {
  inc_mat <- get.edgelist(graph, names = F)
  membership_head <- membership[inc_mat[, 1]]
  membership_tail <- membership[inc_mat[, 2]]
  edge_status <- rep("w", ecount(graph))
  edge_status[membership_head != membership_tail] <- "b"
  return(edge_status)
}

# function to split an existing cluster given MST

splitCluster <- function(mstgraph, k, membership) {
  tcluster <- table(membership)
  clust.split <- sample.int(k, 1, prob = tcluster - 1, replace = TRUE)
  edge_cutted <- sample.int(tcluster[clust.split] - 1, 1)

  mst_subgraph <- igraph::induced_subgraph(mstgraph, membership == clust.split)
  mst_subgraph <- delete.edges(mst_subgraph, edge_cutted)
  connect_comp <- components(mst_subgraph)
  cluster_new <- connect_comp$membership
  vid_new <- (V(mst_subgraph)$vid)[cluster_new == 2] # vid for vertices belonging to new cluster
  membership[vid_new] <- k + 1
  return(list(
    membership = membership, vid_new = vid_new,
    cluster_old = clust.split,
    cluster_new = k + 1
  ))
}



# function to merge two existing clusters

mergeCluster <- function(mstgraph, edge_status, membership) {
  # candidate edges for merging
  ecand <- E(mstgraph)[edge_status == "b"]
  edge_merge <- ecand[sample.int(length(ecand), 1)]
  # update cluster information
  # endpoints of edge_merge, note v1$vid > v2$vid
  v1 <- head_of(mstgraph, edge_merge)
  v2 <- tail_of(mstgraph, edge_merge)
  # clusters that v1, v2 belonging to

  c1 <- membership[v1]
  c2 <- membership[v2]

  ### merge the cluster with a larger label into the one with a smaller label
  if (c1 < c2) {
    c_rm <- c2
    c_new <- c1
  } else {
    c_rm <- c1
    c_new <- c2
  }

  idx_rm <- (membership == c_rm)

  # vid of vertices in c_rm
  vid_old <- (V(mstgraph))[idx_rm]

  # now drop c_rm
  membership[idx_rm] <- c_new
  membership[membership > c_rm] <- membership[membership > c_rm] - 1

  # return the membership, the indices of the merged vertices
  return(list(membership = membership, vid_old = vid_old, cluster_rm = c_rm, cluster_new = c_new))
}


# function to propose a new MST
proposeMST <- function(graph0, edge_status_G) {
  nedge <- length(edge_status_G)
  nb <- sum(edge_status_G == "b")
  nw <- nedge - nb
  weight <- numeric(nedge)
  weight[edge_status_G == "w"] <- runif(nw)
  weight[edge_status_G == "b"] <- runif(nb, 100, 200)
  mstgraph <- mst(graph0, weights = weight)
  return(mstgraph)
}

ord.mat <- function(M, decr = F, cols = NULL) {
  if (is.null(cols)) {
    cols <- 1:ncol(M)
  }
  out <- do.call("order", as.data.frame(M[, cols]))
  if (decr) {
    out <- rev(out)
  }
  return(M[out, ])
}
