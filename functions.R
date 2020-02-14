

# create toy example data
get_dummy_data <- function() {
  layer_suitability_t1 <- raster(
    nrows = 7,
    ncols = 1,
    vals = c(0, 1, 0, 0, 0, 1, 1),
    res = 1,
    xmn = -0.5,
    xmx = 0.5,
    ymn = -4,
    ymx = 3)
  layer_suitability_t2 <- raster(
    nrows = 7,
    ncols = 1,
    vals = c(1, 1, 1, 0, 0, 0, 1),
    res = 1,
    xmn = -0.5,
    xmx = 0.5,
    ymn = -4,
    ymx = 3)
  layer_cost <- raster(
    nrows = 7,
    ncols = 1,
    vals = c(1, 1, 1, NA, NA, 1, 1),
    res = 1,
    xmn = -0.5,
    xmx = 0.5,
    ymn = -4,
    ymx = 3)
  stack_suitability <- stack(
    layer_suitability_t1,
    layer_suitability_t2)
  names(stack_suitability) <- c('T0', 'T1')
  return(stack_suitability, layer_cost)
}


# TODO: understand this function
accCost2 <- function(x, fromCoords) {
  fromCells <- cellFromXY(x, fromCoords)
  tr <- transitionMatrix(x)
  tr <- rbind(tr, rep(0, nrow(tr)))
  tr <- cbind(tr, rep(0, nrow(tr)))
  startNode <- nrow(tr)
  adjP <- cbind(rep(startNode, times = length(fromCells)), fromCells)
  tr[adjP] <- Inf
  adjacencyGraph <- graph.adjacency(tr, mode = 'directed', weighted = TRUE)
  E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight
  return(shortest.paths(adjacencyGraph, v = startNode, mode = 'out')[-startNode])
}
