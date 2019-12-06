### Create and visualize graph

library("igraph")

visualizeGraph <- function(result, rho, lK) {
  w <- result[["W"]]
  n <- length(lK)
  listInd <- rep(1:n, times = lK)
  N <- length(listInd)
  listEdges <- t(matrix(c(listInd, rho), 2,N, byrow = TRUE))
  g <- graph_from_edgelist(listEdges, directed = TRUE)
  listCom <- sapply(1:n, function(i) which.max(w[i, ]))
  print("Here is colour")
  V(g)$color <- listCom
  print("Colour capout")
  add_vertices(g, sum(lK == 0))
  plot(g)
}

visualizeSimpleGraph <- function(rho, lK) {
  n <- length(lK)
  listInd <- rep(1:n, times = lK)
  N <- length(listInd)
  listEdges <- t(matrix(c(listInd, rho), 2,N, byrow = TRUE))
  g <- graph_from_edgelist(listEdges, directed = TRUE)
  g <- add_vertices(g, sum(lK == 0))
  plot(g)
}