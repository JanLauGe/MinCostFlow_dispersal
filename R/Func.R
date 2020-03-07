#' helper function to add new constraints to an already existsing
#' list of constraints. Existing list should be length = 3 with elements
#' - lhs = the left side of the optimization equation (coefficients)
#' - dir = the condition operation, i.e. '>=', '==', or '<='
#' - rhs = the constraint value
#' @param constraints An constraintsMatrix dataset
#' @param new_lhs A new left side of the optimization equation (coefficients)
#' @param new_dir A new vector of condition operations, i.e. c('>=', '==', or '<=')
#' @param new_rhs A new vector of constraint values
#' @return Adds constraints to a constraints matrix object
#' @author Laurens Geffert <laurensgeffert@gmail.com>
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>

add_constraints <- function(
  constraints,
  new_lhs,
  new_dir,
  new_rhs
) {

  #' TODO: add checks for names, length, etc.
  constraints[['lhs']] <- rbind(constraints[['lhs']], new_lhs)
  constraints[['dir']] <- c(constraints[['dir']], new_dir)
  constraints[['rhs']] <- c(constraints[['rhs']], new_rhs)
  return(constraints)
}

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
#' @param edges a data frame with the edge list
#' @param total_flow a numeric value for the total flow
#' @return Creates constraint matrix dataset, this is a list with three elements:
#' - lhs = the left side of the optimization equation (coefficients)
#' - dir = the condition operation, i.e. '>=', '==', or '<='
#' - rhs = the constraint value
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom magrittr "%>%"
#' @importFrom magrittr set_colnames
#' @importFrom magrittr set_rownames
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author Laurens Geffert <laurensgeffert@gmail.com>


create_constraints_matrix <- function(edges, total_flow) {
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

#' Creates a list with the distances among nodes
#' @param x Define this
#' @param fromCoords define that
#' @return Creates a list with the distances among nodes
#' @importFrom igraph E "E<-"
#' @importFrom igraph graph.adjacency
#' @importFrom igraph shortest.paths
#' @importFrom gdistance transitionMatrix
#' @importFrom raster cellFromXY
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author Laurens Geffert <laurensgeffert@gmail.com>

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

#' Creates a data frame with the distances among nodes
#' @param layers_habitat A stack with the binary projection of Species Distribution Models where
#' every layer is a time-slice
#' @param layer_cost A raster of cost of each raster cell in the layers_habitat data
#' @param Dist Maximum distance in meters that the species can travel between time-slices
#' @return Creates a data frame with the distances among nodes keeping only the
#' edges that have a distance shorter than Dist
#' @importFrom gdistance geoCorrection
#' @importFrom gdistance transition
#' @importFrom purrr map_dfr
#' @importFrom raster ncell
#' @importFrom raster values "values<-"
#' @importFrom raster xyFromCell
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author Laurens Geffert <laurensgeffert@gmail.com>

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

#' Solves the linear version of Network Flow
#' @param layers_habitat A stack with the binary projection of Species Distribution Models where
#' every layer is a time-slice
#' @param layer_cost A raster of cost of each raster cell in the layers_habitat data
#' @param Dist Maximum distance in meters that the species can travel between time-slices
#' @param nchains The minimum number of cells needed for the target species on each time-slice
#' @return A list with a rasterstack with the solution time-slice by time-slice and the global
#' solution
#' @examples
#' # Example 1 based on the Phillips problem
#' #Load the Species distribution raster
#' data("Phillips")
#'
#' # see the dataset
#' library(raster)
#' plot(Phillips)
#'
#' # Load the cost layer
#'
#' data("PhilCost")
#'
#' # See the dataset
#'
#' plot(PhilCost)
#'
#' # Solve the linear problem with nchains = 2 and a dispersal of 111 km per time-slice
#'
#' Solution <- solution_linear(layers_habitat = Phillips,
#' layer_cost = PhilCost, Dist = 111000, nchains = 2)
#'
#' # plot the solution for each time-slice
#'
#' plot(Solution$FinalStack)
#'
#' # plot the final solution
#'
#' plot(Solution$FinalSolution)
#'
#' #' # Example 2 with hypothetical species
#' #Load the Species distribution raster
#' data("BinSpp")
#'
#' # see the dataset
#'
#' plot(BinSpp[[1]])
#'
#' # Load the cost layer
#'
#' data("Cost")
#'
#' # See the dataset
#'
#' plot(Cost)
#'
#' # Solve the linear problem with nchains = 2 and a dispersal of 111 km per time-slice
#'
#' Solution <- solution_linear(layers_habitat = BinSpp[[1]],
#' layer_cost = Cost, Dist = 240000, nchains = 5)
#'
#' # plot the solution for each time-slice
#'
#' plot(Solution$FinalStack)
#'
#' # plot the final solution
#'
#' plot(Solution$FinalSolution)
#'
#' @importFrom dplyr arrange
#' @importFrom dplyr bind_cols
#' @importFrom dplyr bind_rows
#' @importFrom dplyr group_by
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr rename
#' @importFrom dplyr row_number
#' @importFrom dplyr select
#' @importFrom dplyr summarise
#' @importFrom dplyr tibble
#' @importFrom dplyr transmute
#' @importFrom magrittr "%>%"
#' @importFrom lpSolve lp
#' @importFrom raster nlayers
#' @importFrom raster values "values<-"
#' @importFrom tidyr expand_grid
#' @importFrom purrr map_dfr
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author Laurens Geffert <laurensgeffert@gmail.com>
#' @export

solution_linear <- function(layers_habitat, layer_cost, Dist, nchains){
  layers_habitat <- layers_habitat * !is.na(layer_cost)

  df_habitat <- map_dfr(
    layers_habitat = layers_habitat,
    .x = 1:nlayers(layers_habitat),
    .f = ~layers_habitat[[.x]] %>%
      values() %>%
      tibble(habitat = .) %>%
      transmute(
        cell_id = row_number(),
        time = .x,
        habitat = habitat)) %>%
    # give the same cell in different time slices different IDs
    mutate(raster_id = cell_id, cell_id = cell_id + max(cell_id) * time)

  df_IDs <- df_habitat %>% rename(node_to = cell_id) %>% dplyr::select(node_to, raster_id)

  df_cost <- layer_cost %>%
    values() %>%
    tibble(cost = .) %>%
    mutate(cell_id = row_number())
  df_cost <- map_dfr(
    df_cost = df_cost,
    .x = 1:nlayers(layers_habitat),
    .f = ~df_cost %>%
      transmute(
        cell_id = cell_id + max(cell_id) * .x,
        edge_cost = cost))

  connections <- edge_distance_limit(layer_cost = layer_cost, layers_habitat = layers_habitat, Dist = Dist)

  edges_timesteps <- map_dfr(
    df_habitat = df_habitat,
    .x = 2:nlayers(layers_habitat),
    .f = ~expand_grid(
      node_from = df_habitat %>%
        filter(time == .x - 1 & habitat == 1) %>%
        pull(cell_id),
      node_to = df_habitat %>%
        filter(time == .x & habitat == 1) %>%
        pull(cell_id),
      timeslice_to = .x)) %>% mutate(
        edge_id = row_number(),
        edge_capacity = 1,
        edge_cost = NA)

  for(i in 1:nrow(edges_timesteps)){
    #Check if the raster ID is equal in node_to and node_from
    if(df_IDs$raster_id[df_IDs$node_to == edges_timesteps[i,]$node_to] == df_IDs$raster_id[df_IDs$node_to == edges_timesteps[i,]$node_from]){
      ## If that is the case just add the cost once
      edges_timesteps$edge_cost[i] <- df_cost$edge_cost[df_cost$cell_id == edges_timesteps[i,]$node_from]
    }
    #Check if the raster ID is Different in node_to and node_from
    if(df_IDs$raster_id[df_IDs$node_to == edges_timesteps[i,]$node_to] != df_IDs$raster_id[df_IDs$node_to == edges_timesteps[i,]$node_from]){
      ## If that is the case just add the cost twice
      edges_timesteps$edge_cost[i] <- df_cost$edge_cost[df_cost$cell_id == edges_timesteps[i,]$node_from] + df_cost$edge_cost[df_cost$cell_id == edges_timesteps[i,]$node_to]
    }
  }

  edges_timesteps <- left_join(edges_timesteps, connections) %>% dplyr::filter(!is.na(dist)) %>% dplyr::select(-dist)

  edges_source <- tibble(
    # lowest ids available
    # FIXME: temp workaround
    edge_id = min(edges_timesteps[['edge_id']]),
    node_from = 0,
    # find all timeslice 1 nodes
    node_to = edges_timesteps %>%
      filter(timeslice_to == 2) %>%
      pull(node_from) %>%
      unique(),
    timeslice_to = 1,
    edge_capacity = nchains,
    edge_cost = 1) %>%
    # unique ids!
    mutate(edge_id = edge_id - row_number())

  edges_target <- tibble(
    # FIXME: temp workaround
    edge_id = max(edges_timesteps[['edge_id']]),
    node_from = edges_timesteps %>%
      filter(timeslice_to == max(edges_timesteps[['timeslice_to']])) %>%
      pull(node_to) %>%
      unique(),
    node_to = max(edges_timesteps[['node_to']]) + 1,
    timeslice_to = max(edges_timesteps[['timeslice_to']]) + 1,
    edge_capacity = nchains,
    edge_cost = 1) %>%
    # unique ids!
    mutate(edge_id = edge_id + row_number())

  edges_formatted <- edges_timesteps %>%
    bind_rows(edges_source, edges_target) %>%
    arrange(edge_id)

  constraints_matrix <- create_constraints_matrix(
    edges = edges_formatted,
    total_flow = nchains)

  solution <- lp(
    direction = 'min',
    objective.in = edges_formatted[['edge_cost']],
    const.mat = constraints_matrix[['lhs']],
    const.dir = constraints_matrix[['dir']],
    const.rhs = constraints_matrix[['rhs']])

  edges_solved <- bind_cols(edges_formatted, solution = solution[['solution']]) %>% left_join(df_IDs) %>% dplyr::filter(solution != 0)


  FinalStack <- list()


  for(i in 1:nlayers(layers_habitat)){
    Temp<- layer_cost
    values(Temp) <- 0
    values(Temp)[edges_solved %>% dplyr::filter(timeslice_to == (i)) %>% pull(raster_id) %>%  unique()] <- edges_solved %>% dplyr::filter(timeslice_to == (i)) %>% group_by(raster_id) %>% summarise(flow = sum(solution))  %>% pull(flow)
    FinalStack[[i]] <- Temp
  }

  FinalStack <- do.call("stack", FinalStack)
  names(FinalStack) <- paste0("T", 1:nlayers(FinalStack))

  FinalSolution <- max(FinalStack)

  return(list(FinalStack = FinalStack, FinalSolution = FinalSolution))
}

#' Solves the linear version of Network Flow using the old algorithm
#' @param layers_habitat A stack with the binary projection of Species Distribution Models where
#' every layer is a time-slice
#' @param layer_cost A raster of cost of each raster cell in the layers_habitat data
#' @param Dist Maximum distance in meters that the species can travel between time-slices
#' @param nchains The minimum number of cells needed for the target species on each time-slice
#' @return A list with a rasterstack with the solution time-slice by time-slice and the global
#' solution
#' @examples
#' # Example 1 based on the Phillips problem
#' #Load the Species distribution raster
#' data("Phillips")
#'
#' # see the dataset
#' library(raster)
#' plot(Phillips)
#'
#' # Load the cost layer
#'
#' data("PhilCost")
#'
#' # See the dataset
#'
#' plot(PhilCost)
#'
#' # Solve the linear problem with nchains = 2 and a dispersal of 111 km per time-slice
#'
#' Solution <- solution_linear_old(layers_habitat = Phillips,
#' layer_cost = PhilCost, Dist = 111000, nchains = 2)
#'
#' # plot the solution for each time-slice
#'
#' plot(Solution$FinalStack)
#'
#' # plot the final solution
#'
#' plot(Solution$FinalSolution)
#'
#' #' # Example 2 with hypothetical species
#' #Load the Species distribution raster
#' data("BinSpp")
#'
#' # see the dataset
#'
#' plot(BinSpp[[1]])
#'
#' # Load the cost layer
#'
#' data("Cost")
#'
#' # See the dataset
#'
#' plot(Cost)
#'
#' # Solve the linear problem with nchains = 2 and a dispersal of 111 km per time-slice
#'
#' Solution <- solution_linear_old(layers_habitat = BinSpp[[1]],
#' layer_cost = Cost, Dist = 240000, nchains = 5)
#'
#' # plot the solution for each time-slice
#'
#' plot(Solution$FinalStack)
#'
#' # plot the final solution
#'
#' plot(Solution$FinalSolution)
#'
#' @importFrom dplyr arrange
#' @importFrom dplyr bind_cols
#' @importFrom dplyr bind_rows
#' @importFrom dplyr group_by
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr rename
#' @importFrom dplyr row_number
#' @importFrom dplyr select
#' @importFrom dplyr summarise
#' @importFrom dplyr tibble
#' @importFrom dplyr transmute
#' @importFrom magrittr "%>%"
#' @importFrom lpSolve lp
#' @importFrom raster nlayers
#' @importFrom raster values "values<-"
#' @importFrom tidyr expand_grid
#' @importFrom purrr map_dfr
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author Laurens Geffert <laurensgeffert@gmail.com>
#' @export

solution_linear_old <- function(layers_habitat, layer_cost, Dist, nchains){
  layers_habitat <- layers_habitat * !is.na(layer_cost)

  df_habitat <- map_dfr(
    layers_habitat = layers_habitat,
    .x = 1:nlayers(layers_habitat),
    .f = ~layers_habitat[[.x]] %>%
      values() %>%
      tibble(habitat = .) %>%
      transmute(
        cell_id = row_number(),
        time = .x,
        habitat = habitat)) %>%
    # give the same cell in different time slices different IDs
    mutate(raster_id = cell_id, cell_id = cell_id + max(cell_id) * time)

  df_IDs <- df_habitat %>% rename(node_to = cell_id) %>% dplyr::select(node_to, raster_id)

  df_cost <- layer_cost %>%
    values() %>%
    tibble(cost = .) %>%
    mutate(cell_id = row_number())
  df_cost <- map_dfr(
    df_cost = df_cost,
    .x = 1:nlayers(layers_habitat),
    .f = ~df_cost %>%
      transmute(
        cell_id = cell_id + max(cell_id) * .x,
        edge_cost = cost))

  connections <- edge_distance_limit(layer_cost = layer_cost, layers_habitat = layers_habitat, Dist = Dist)

  edges_timesteps <- map_dfr(
    df_habitat = df_habitat,
    .x = 2:nlayers(layers_habitat),
    .f = ~expand_grid(
      node_from = df_habitat %>%
        filter(time == .x - 1 & habitat == 1) %>%
        pull(cell_id),
      node_to = df_habitat %>%
        filter(time == .x & habitat == 1) %>%
        pull(cell_id),
      timeslice_to = .x)) %>% mutate(
        edge_id = row_number(),
        edge_capacity = 1,
        edge_cost = NA)

  for(i in 1:nrow(edges_timesteps)){
    edges_timesteps$edge_cost[i] <- df_cost$edge_cost[df_cost$cell_id == edges_timesteps[i,]$node_from] + df_cost$edge_cost[df_cost$cell_id == edges_timesteps[i,]$node_to]
  }

  edges_timesteps <- left_join(edges_timesteps, connections) %>% dplyr::filter(!is.na(dist)) %>% dplyr::select(-dist)

  edges_source <- tibble(
    # lowest ids available
    # FIXME: temp workaround
    edge_id = min(edges_timesteps[['edge_id']]),
    node_from = 0,
    # find all timeslice 1 nodes
    node_to = edges_timesteps %>%
      filter(timeslice_to == 2) %>%
      pull(node_from) %>%
      unique(),
    timeslice_to = 1,
    edge_capacity = nchains,
    edge_cost = 1) %>%
    # unique ids!
    mutate(edge_id = edge_id - row_number())

  edges_target <- tibble(
    # FIXME: temp workaround
    edge_id = max(edges_timesteps[['edge_id']]),
    node_from = edges_timesteps %>%
      filter(timeslice_to == max(edges_timesteps[['timeslice_to']])) %>%
      pull(node_to) %>%
      unique(),
    node_to = max(edges_timesteps[['node_to']]) + 1,
    timeslice_to = max(edges_timesteps[['timeslice_to']]) + 1,
    edge_capacity = nchains,
    edge_cost = 1) %>%
    # unique ids!
    mutate(edge_id = edge_id + row_number())

  edges_formatted <- edges_timesteps %>%
    bind_rows(edges_source, edges_target) %>%
    arrange(edge_id)

  constraints_matrix <- create_constraints_matrix(
    edges = edges_formatted,
    total_flow = nchains)

  solution <- lp(
    direction = 'min',
    objective.in = edges_formatted[['edge_cost']],
    const.mat = constraints_matrix[['lhs']],
    const.dir = constraints_matrix[['dir']],
    const.rhs = constraints_matrix[['rhs']])

  edges_solved <- bind_cols(edges_formatted, solution = solution[['solution']]) %>% left_join(df_IDs) %>% dplyr::filter(solution != 0)


  FinalStack <- list()


  for(i in 1:nlayers(layers_habitat)){
    Temp<- layer_cost
    values(Temp) <- 0
    values(Temp)[edges_solved %>% dplyr::filter(timeslice_to == (i)) %>% pull(raster_id) %>%  unique()] <- edges_solved %>% dplyr::filter(timeslice_to == (i)) %>% group_by(raster_id) %>% summarise(flow = sum(solution))  %>% pull(flow)
    FinalStack[[i]] <- Temp
  }

  FinalStack <- do.call("stack", FinalStack)
  names(FinalStack) <- paste0("T", 1:nlayers(FinalStack))

  FinalSolution <- max(FinalStack)

  return(list(FinalStack = FinalStack, FinalSolution = FinalSolution))
}


#' Solves the quadratic version of Network Flow
#' @param layers_habitat A stack with the binary projection of Species Distribution Models where
#' every layer is a time-slice
#' @param layer_cost A raster of cost of each raster cell in the layers_habitat data
#' @param Dist Maximum distance in meters that the species can travel between time-slices
#' @param nchains The minimum number of cells needed for the target species on each time-slice
#' @return A list with a rasterstack with the solution time-slice by time-slice and the global
#' solution
#' @examples
#'  # Example 1 based on the Phillips problem
#'  #Load the Species distribution raster
#' data("Phillips")
#'
#' # see the dataset
#' library(raster)
#' plot(Phillips)
#'
#' # Load the cost layer
#'
#' data("PhilCost")
#'
#' # See the dataset
#'
#' plot(PhilCost)
#'
#' # Solve the linear problem with nchains = 2 and a dispersal of 111 km per time-slice
#'
#' Solution <- solution_quadratic(layers_habitat = Phillips,
#' layer_cost = PhilCost, Dist = 111000, nchains = 2)
#'
#' # plot the solution for each time-slice
#'
#' plot(Solution$FinalStack)
#'
#' # plot the final solution
#'
#' plot(Solution$FinalSolution)
#'
#' #' # Example 2 with hypothetical species
#' #Load the Species distribution raster
#' data("BinSpp")
#'
#' # see the dataset
#'
#' plot(BinSpp[[1]])
#'
#' # Load the cost layer
#'
#' data("Cost")
#'
#' # See the dataset
#'
#' plot(Cost)
#'
#' # Solve the linear problem with nchains = 2 and a dispersal of 111 km per time-slice
#'
#' Solution <- solution_quadratic(layers_habitat = BinSpp[[1]],
#' layer_cost = Cost, Dist = 240000, nchains = 5)
#'
#' # plot the solution for each time-slice
#'
#' plot(Solution$FinalStack)
#'
#' # plot the final solution
#'
#' plot(Solution$FinalSolution)
#'
#' @importFrom dplyr arrange
#' @importFrom dplyr bind_cols
#' @importFrom dplyr bind_rows
#' @importFrom dplyr group_by
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr rename
#' @importFrom dplyr row_number
#' @importFrom dplyr select
#' @importFrom dplyr summarise
#' @importFrom dplyr tibble
#' @importFrom dplyr transmute
#' @importFrom magrittr "%>%"
#' @importFrom optiSolve cop
#' @importFrom optiSolve lbcon
#' @importFrom optiSolve lincon
#' @importFrom optiSolve quadfun
#' @importFrom optiSolve solvecop
#' @importFrom optiSolve ubcon
#' @importFrom raster nlayers
#' @importFrom raster values "values<-"
#' @importFrom tidyr expand_grid
#' @importFrom purrr map_dfr
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author Laurens Geffert <laurensgeffert@gmail.com>
#' @export

solution_quadratic <- function(layers_habitat, layer_cost, Dist, nchains){
  layers_habitat <- layers_habitat * !is.na(layer_cost)

  df_habitat <- map_dfr(
    layers_habitat = layers_habitat,
    .x = 1:nlayers(layers_habitat),
    .f = ~layers_habitat[[.x]] %>%
      values() %>%
      tibble(habitat = .) %>%
      transmute(
        cell_id = row_number(),
        time = .x,
        habitat = habitat)) %>%
    # give the same cell in different time slices different IDs
    mutate(raster_id = cell_id, cell_id = cell_id + max(cell_id) * time)

  df_IDs <- df_habitat %>% rename(node_to = cell_id) %>% dplyr::select(node_to, raster_id)

  df_cost <- layer_cost %>%
    values() %>%
    tibble(cost = .) %>%
    mutate(cell_id = row_number())
  df_cost <- map_dfr(
    df_cost = df_cost,
    .x = 1:nlayers(layers_habitat),
    .f = ~df_cost %>%
      transmute(
        cell_id = cell_id + max(cell_id) * .x,
        edge_cost = cost))

  connections <- edge_distance_limit(layer_cost = layer_cost, layers_habitat = layers_habitat, Dist = Dist)

  edges_timesteps <- map_dfr(
    df_habitat = df_habitat,
    .x = 2:nlayers(layers_habitat),
    .f = ~expand_grid(
      node_from = df_habitat %>%
        filter(time == .x - 1 & habitat == 1) %>%
        pull(cell_id),
      node_to = df_habitat %>%
        filter(time == .x & habitat == 1) %>%
        pull(cell_id),
      timeslice_to = .x)) %>% mutate(
        edge_id = row_number(),
        edge_capacity = 1,
        edge_cost = NA)

  for(i in 1:nrow(edges_timesteps)){
    #Check if the raster ID is equal in node_to and node_from
    if(df_IDs$raster_id[df_IDs$node_to == edges_timesteps[i,]$node_to] == df_IDs$raster_id[df_IDs$node_to == edges_timesteps[i,]$node_from]){
      ## If that is the case just add the cost once
      edges_timesteps$edge_cost[i] <- df_cost$edge_cost[df_cost$cell_id == edges_timesteps[i,]$node_from]
    }
    #Check if the raster ID is Different in node_to and node_from
    if(df_IDs$raster_id[df_IDs$node_to == edges_timesteps[i,]$node_to] != df_IDs$raster_id[df_IDs$node_to == edges_timesteps[i,]$node_from]){
      ## If that is the case just add the cost twice
      edges_timesteps$edge_cost[i] <- df_cost$edge_cost[df_cost$cell_id == edges_timesteps[i,]$node_from] + df_cost$edge_cost[df_cost$cell_id == edges_timesteps[i,]$node_to]
    }
  }

  edges_timesteps <- left_join(edges_timesteps, connections) %>% dplyr::filter(!is.na(dist)) %>% dplyr::select(-dist)

  edges_source <- tibble(
    # lowest ids available
    # FIXME: temp workaround
    edge_id = min(edges_timesteps[['edge_id']]),
    node_from = 0,
    # find all timeslice 1 nodes
    node_to = edges_timesteps %>%
      filter(timeslice_to == 2) %>%
      pull(node_from) %>%
      unique(),
    timeslice_to = 1,
    edge_capacity = nchains,
    edge_cost = 1) %>%
    # unique ids!
    mutate(edge_id = edge_id - row_number())

  edges_target <- tibble(
    # FIXME: temp workaround
    edge_id = max(edges_timesteps[['edge_id']]),
    node_from = edges_timesteps %>%
      filter(timeslice_to == max(edges_timesteps[['timeslice_to']])) %>%
      pull(node_to) %>%
      unique(),
    node_to = max(edges_timesteps[['node_to']]) + 1,
    timeslice_to = max(edges_timesteps[['timeslice_to']]) + 1,
    edge_capacity = nchains,
    edge_cost = 1) %>%
    # unique ids!
    mutate(edge_id = edge_id + row_number())

  edges_formatted <- edges_timesteps %>%
    bind_rows(edges_source, edges_target) %>%
    arrange(edge_id)

  constraints_matrix <- create_constraints_matrix(
    edges = edges_formatted,
    total_flow = nchains)

  DmatTest <- matrix(0, ncol = length(edges_formatted$edge_cost), nrow = length(edges_formatted$edge_cost))


  for(i in 1:length(edges_formatted$edge_cost)){
    DmatTest[i,i] <- (edges_formatted$edge_cost[i])
  }

  mycop <- cop(f  = quadfun(Q = DmatTest),
               lb = lbcon(0),
               ub = ubcon(NA),
               lc = lincon(A=constraints_matrix$lhs, dir=constraints_matrix$dir, val=constraints_matrix$rhs, name = 1:nrow(constraints_matrix$lhs)))
  suppressWarnings({
    res <- solvecop(mycop, solver="cccp", quiet=FALSE, trace=FALSE)
  })

  edges_solved <- bind_cols(edges_formatted, solution = res$x) %>% left_join(df_IDs) %>% dplyr::filter(solution != 0) %>% group_by(node_to, raster_id, timeslice_to) %>% summarise(solution = sum(solution)) %>% mutate(solution = ifelse(solution < 0, 0, solution))


  FinalStack <- list()


  for(i in 1:nlayers(layers_habitat)){
    Temp<- layer_cost
    values(Temp) <- 0
    values(Temp)[edges_solved %>% dplyr::filter(timeslice_to == (i)) %>% pull(raster_id) %>%  unique()] <- edges_solved %>% dplyr::filter(timeslice_to == (i)) %>% group_by(raster_id) %>% summarise(flow = sum(solution))  %>% pull(flow)
    FinalStack[[i]] <- Temp
  }

  FinalStack <- do.call("stack", FinalStack)
  names(FinalStack) <- paste0("T", 1:nlayers(FinalStack))

  FinalSolution <- max(FinalStack)

  return(list(FinalStack = FinalStack, FinalSolution = FinalSolution, edges_solved = edges_solved))
}
