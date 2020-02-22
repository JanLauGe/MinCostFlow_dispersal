

get_dummy_data <- function() {
  #' helper function that creates example raster data
  #' consisting of habitat suitability layers for two time steps
  #' as well as a cost layer
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


data_phillips_problem <- function() {
  #' helper function to re-create phillips problem in edge list format
  edges_phillips <- data.frame(
    edge_id = seq(1, 14, 1),
    node_from = c(1, 1, 1, 1, 2, 2, 2, 3, 4, 5, 6, 7, 8, 9),
    node_to = c(2, 3, 4, 5, 6, 7, 8, 9, 9, 9, 10, 10, 10, 10),
    edge_capacity = rep(1, 14),
    edge_cost = rep(1, 14))
  return(edges_phillips)
}


add_constraints <- function(
  constraints,
  new_lhs,
  new_dir,
  new_rhs
) {
  #' helper function to add new constraints to an already existsing
  #' list of constraints. Existing list should be length = 3 with elements
  #' - lhs = the left side of the optimization equation (coefficients)
  #' - dir = the condition operation, i.e. '>=', '==', or '<='
  #' - rhs = the constraint value
  #' TODO: add checks for names, length, etc.
  constraints[['lhs']] <- rbind(constraints[['lhs']], new_lhs)
  constraints[['dir']] <- c(constraints[['dir']], new_dir)
  constraints[['rhs']] <- c(constraints[['rhs']], new_rhs)
  return(constraints)
}


create_constraints_matrix <- function(edges, total_flow) {
  #' creates a constraint matrix for solving a minimum cost flow
  #' problem given a set of edges as data frame in the following format:
  #' - edge_id (primary key)
  #' - node_from (foreign key)
  #' - node_to (foreign key)
  #' - edge_capacity (int)
  #' - edge_cost (int)
  #' Assumptions:
  #' - source node is the node with the lowest ID
  #' - target node is the node with the highest ID

  # extract information on edges
  ids_edges <- edges[['edge_id']]
  n_edges <- length(ids_edges)
  # extract information on nodes
  # Source node is assumed to be the one with the minimum ID number!
  # Sink node is assumed to be the one with the maximum ID number!
  node_id_source <- min(edges[['node_from']])
  node_id_target <- max(edges[['node_to']])
  ids_nodes <- c(edges[['node_from']], edges[['node_to']]) %>% unique()
  n_nodes <- length(ids_nodes)
  # initialize empty constraints matrix
  constraints <- list(
    lhs = NA,
    dir = NA,
    rhs = NA)

  # build edge capacity constraints --------------------------------------------
  #' Flow through each edge should not be larger than capacity.
  #' We create one constraint for each edge. All coefficients zero
  #' except the ones of the edge in question as one, with a constraint
  #' that the result is smaller than or equal to capacity of that edge.
  #' TODO: We may be able to take this out if we are using node capacity constraints
  constraints[['lhs']] <- edges %>%
    pull(edge_id) %>%
    length() %>%
    diag() %>%
    set_colnames(edges[['edge_id']]) %>%
    set_rownames(edges[['edge_id']])
  # should be smaller than or equal to
  constraints[['dir']] <- rep('<=', times = nrow(edges))
  # than capacity
  constraints[['rhs']] <- edges[['edge_capacity']]

  # build node residual constraints ------------------------------------------------
  #' No node except for the source and target node should retain any flow.
  #' Therefore, the sum of all inputs and outputs for these nodes should be
  #' zero. We can force this by going through all nodes and identifying their
  #' edges. The inbound edges are assigned a coefficient of 1 and the outbound
  #' edges a coefficients of -1. Any viable solution should achieve a result
  #' equal to zero, i.e. inputs equal to outputs for all nodes (except S and T).
  constraint_node_residuals <- matrix(
    data = 0,
    nrow = n_nodes,
    ncol = n_edges,
    dimnames = list(ids_nodes, ids_edges))
  # loop over nodes
  for (i in ids_nodes) {
    # input edges
    edges_in <- edges %>%
      filter(node_to == i) %>%
      pull(edge_id)
    # output edges
    edges_out <- edges %>%
      filter(node_from == i) %>%
      pull(edge_id)
    # set input coefficients to 1
    constraint_node_residuals[
      rownames(constraint_node_residuals) == i,
      colnames(constraint_node_residuals) %in% edges_in] <- 1
    # set output coefficients to -1
    constraint_node_residuals[
      rownames(constraint_node_residuals) == i,
      colnames(constraint_node_residuals) %in% edges_out] <- -1
  }
  # exclude source and target edges.
  # keep node flow values for separate step below
  node_flow_source <- constraint_node_residuals[rownames(constraint_node_residuals) == node_id_source,]
  node_flow_target <- constraint_node_residuals[rownames(constraint_node_residuals) == node_id_target,]
  # exclude them from node flow here
  constraint_node_residuals <- constraint_node_residuals[!rownames(constraint_node_residuals) %in% c(node_id_source, node_id_target),]
  # add node residual constraints to the constraints list
  constraints <- add_constraints(
    constraints = constraints,
    new_lhs = constraint_node_residuals,
    new_dir = rep('==', times = nrow(constraint_node_residuals)),
    new_rhs = rep(0, times = nrow(constraint_node_residuals)))

  # build node capacity constraints --------------------------------------------
  #' In this particular case we want the flow through each node to be smaller
  #' than or equal to 1 so that the number of chains is preserved. We add this
  #' constraint in a very similar way to the node flow constraints
  # take previously computed information on incoming nodes
  contraint_node_capacity <- constraint_node_residuals
  # discard information on outgoing nodes
  contraint_node_capacity[contraint_node_capacity < 1] <- 0
  # add node capacity constraints to the constraints list
  constraints <- add_constraints(
    constraints = constraints,
    new_lhs = contraint_node_capacity,
    new_dir = rep('<=', times = nrow(contraint_node_capacity)),
    new_rhs = rep(1, times = nrow(contraint_node_capacity)))

  # Build initialisation constraints -------------------------------------------
  #' For the source and the target node, we want all outbound nodes and
  #' all inbound nodes to be equal to the sum of flow through the network
  #' respectively
  constraints <- add_constraints(
    constraints = constraints,
    new_lhs = rbind(source = node_flow_source, target = node_flow_target),
    new_dir = rep('==', times = 2),
    new_rhs = c(total_flow * -1, total_flow))

  return(constraints)
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

edge_distance_limit <- function(layer_cost, layers_habitat, Dist){
  Masklayer <- layer_cost
  values(Masklayer) <- ifelse(is.na(values(Masklayer)), NA, 1)

  Stack <- layers_habitat * Masklayer

  Raster <- sum(Stack)
  Raster[values(Raster) > 0] = 1
  Raster[values(Raster) == 0] = NA

  h16 <- transition(Raster, transitionFunction=function(x){1},16,symm=FALSE)
  h16 <- geoCorrection(h16, scl=FALSE)
  ID <-c(1:ncell(Raster))[!is.na(values(Raster))]

  B <- xyFromCell(Raster, cell = ID)
  connections <- list()
  #For each pair of cells in B
  for (i in 1:nrow(B)){
    #Create a temporal raster for each row with the distance from cell xy to all other cells
    temp <- accCost2(h16,B[i,])
    index <- which(temp < Dist)
    connections[[i]] <- cbind(ID[i], index, temp[index])
  }


  connections <- do.call("rbind", connections)
  connections <- as.data.frame(connections)

  colnames(connections) <- c("from", "to", "dist")

  connections <- map_dfr(
    connections = connections,
    .x = 1:nlayers(layers_habitat),
    .f = ~connections %>%
      mutate(
        node_from = from + ncell(layer_cost) * .x,
        node_to= ncell(layer_cost) + (to + ncell(layer_cost) * .x)))

  connections <- connections %>% dplyr::select(-from, -to)
  return(connections)
}
