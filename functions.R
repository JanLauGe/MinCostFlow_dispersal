

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
  layers_habitat <- stack(
    layer_suitability_t1,
    layer_suitability_t2)
  names(layers_habitat) <- c('T0', 'T1')
  return(list(layers_habitat, layer_cost))
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



create_constraints_matrix <- function(edges, total_flow) {
  
  # Edge IDs to be used as names
  names_edges <- edges$ID
  # Number of edges
  numberof_edges <- length(names_edges)
  
  # Node IDs to be used as names
  names_nodes <- c(edges$Cell_From, edges$Cell_To) %>% unique
  # Number of nodes
  numberof_nodes <- length(names_nodes)
  
  # Times id
  
  name_times <- c(edges$Time_From, edges$Time_to) %>% unique
  # Number of nodes
  numberof_times <- length(name_times)
  
  # Build constraints matrix
  constraints <- list(
    lhs = NA,
    dir = NA,
    rhs = NA)
  
  #' Build capacity constraints ------------------------------------------------
  #' Flow through each edge should not be larger than capacity.
  #' We create one constraint for each edge. All coefficients zero
  #' except the ones of the edge in question as one, with a constraint
  #' that the result is smaller than or equal to capacity of that edge.
  
  # Flow through individual edges
  constraints$lhs <- edges$ID %>%
    length %>%
    diag %>%
    magrittr::set_colnames(edges$ID) %>%
    magrittr::set_rownames(edges$ID)
  # should be smaller than or equal to
  constraints$dir <- rep('<=', times = nrow(edges))
  # than capacity
  constraints$rhs <- edges$Capacity
  
  
  #' Build node flow constraints -----------------------------------------------
  #' For each node, find all edges that go to that node
  #' and all edges that go from that node. The sum of all inputs
  #' and all outputs should be zero. So we set inbound edge coefficients as 1
  #' and outbound coefficients as -1. In any viable solution the result should
  #' be equal to zero.
  
  nodeflow <- matrix(0,
                     nrow = numberof_nodes,
                     ncol = numberof_edges,
                     dimnames = list(names_nodes, names_edges))
  
  for (i in names_nodes) {
    # input arcs
    edges_in <- edges %>%
      dplyr::filter(Cell_To == i) %>%
      dplyr::select(ID) %>%
      unlist
    # output arcs
    edges_out <- edges %>%
      dplyr::filter(Cell_From == i) %>%
      dplyr::select(ID) %>%
      unlist
    
    # output arcs
    
    # set input coefficients to 1
    nodeflow[
      rownames(nodeflow) == i,
      colnames(nodeflow) %in% edges_in] <- 1
    
    # set output coefficients to -1
    nodeflow[
      rownames(nodeflow) == i,
      colnames(nodeflow) %in% edges_out] <- -1
  }
  
  # But exclude source and target edges
  # as the zero-sum flow constraint does not apply to these!
  # Source node is assumed to be the one with the minimum ID number
  # Sink node is assumed to be the one with the maximum ID number
  sourcenode_id <- min(edges$Cell_From)
  targetnode_id <- max(edges$Cell_To)
  # Keep node flow values for separate step below
  nodeflow_source <- nodeflow[rownames(nodeflow) == sourcenode_id,]
  nodeflow_target <- nodeflow[rownames(nodeflow) == targetnode_id,]
  # Exclude them from node flow here
  nodeflow <- nodeflow[!rownames(nodeflow) %in% c(sourcenode_id, targetnode_id),]
  
  # Add nodeflow to the constraints list
  constraints$lhs <- rbind(constraints$lhs, nodeflow)
  constraints$dir <- c(constraints$dir, rep('==', times = nrow(nodeflow)))
  constraints$rhs <- c(constraints$rhs, rep(0, times = nrow(nodeflow)))
  
  
  #' Build initialisation constraints ------------------------------------------
  #' For the source and the target node, we want all outbound nodes and
  #' all inbound nodes to be equal to the sum of flow through the network
  #' respectively
  
  # Add initialisation to the constraints list
  constraints$lhs <- rbind(constraints$lhs,
                           source = nodeflow_source,
                           target = nodeflow_target)
  constraints$dir <- c(constraints$dir, rep('==', times = 2))
  # Flow should be negative for source, and positive for target
  constraints$rhs <- c(constraints$rhs, total_flow * -1, total_flow)
  
  ##Generate the time constraints, find all the nodes en each time
  
  Timeflow <- matrix(0,
                     nrow = numberof_times,
                     ncol = numberof_edges,
                     dimnames = list(name_times, names_edges))
  
  for (i in name_times) {
    # input arcs
    edges_in <- edges %>%
      dplyr::filter(Time_to == i) %>%
      dplyr::select(ID) %>%
      unlist
    # output arcs
    edges_out <- edges %>%
      dplyr::filter(Time_From == i) %>%
      dplyr::select(ID) %>%
      unlist
    
    # output arcs
    
    # set input coefficients to 1
    Timeflow[
      rownames(Timeflow) == i,
      colnames(Timeflow) %in% edges_in] <- 1
    
    # set output coefficients to -1
    Timeflow[
      rownames(Timeflow) == i,
      colnames(Timeflow) %in% edges_out] <- -1
  }
  
  # But exclude source and target edges
  # as the zero-sum flow constraint does not apply to these!
  # Source node is assumed to be the one with the minimum ID number
  # Sink node is assumed to be the one with the maximum ID number
  sourcetime_id <- min(edges$Time_From)
  targettime_id <- max(edges$Time_to)
  # Keep node flow values for separate step below
  timeflow_source <- Timeflow[rownames(Timeflow) == sourcetime_id,]
  timeflow_target <- Timeflow[rownames(Timeflow) == targettime_id,]
  # Exclude them from node flow here
  Timeflow <- Timeflow[!rownames(Timeflow) %in% c(sourcetime_id, targettime_id),]
  
  constraints$lhs <- rbind(constraints$lhs, Timeflow)
  constraints$dir <- c(constraints$dir, rep('==', times = nrow(Timeflow)))
  constraints$rhs <- c(constraints$rhs, rep(0, times = nrow(Timeflow)))
  
  ###################
  return(constraints)
}