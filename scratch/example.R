# dependencies
if (!require('pacman')) install.packages('pacman')
if (!require('devtools')) install.packages('devtools')
if (!require('QuadCostAmpl')) devtools::install_github('derek-corcoran-barrios/QuadCostAmpl')
pacman::p_load(devtools, tidyverse, magrittr, gdistance, raster, lpSolve, QuadCostAmpl, tictoc)

source('functions.R')
source('R/funs_min_cost_flow.R')

# example data
data('BinSpp')
data('Cost')

# let's have a look
layers_habitat = BinSpp[[1]]
layer_cost = Cost
Dist = 400000
nchains = 4

plot(layers_habitat)
plot(layer_cost)


NewFunc <- function(
  layers_habitat, # layer_habitat?
  layer_cost, # layer_cost?
  Dist = 400000,
  nchains = 8
) {
  # filter habitat suitability layer for sites with dispersal cost
  layers_habitat <- layers_habitat * !is.na(layer_cost)

  # create dataframe of habitat cells across time slices
  df_habitat <- map_dfr(
    layers_habitat = layers_habitat,
    .x = 1:nlayers(layers_habitat),
    .f = ~layers_habitat[[.x]] %>%
      values() %>%
      tibble(habitat = .) %>%
      transmute(
        id = row_number(),
        time = .x,
        habitat = habitat
      )
  )

  # find cells with at least 1 suitability
  df_habitat <- df_habitat %>%
    group_by(id) %>%
    filter(any(habitat > 0)) %>%
    ungroup()

  # make binary raster
  # FIXME: rename raster to avoid caps and reserved name
  # FIXME: is there a way to do this with the sp package?
  Raster <- sum(layers_habitat)
  Raster[values(Raster) > 0] = 1
  Raster[values(Raster) == 0] = NA

  # calculate the distance between cells across time slices?
  h16 <- Raster %>%
    transition(transitionFunction = function(x) {1}, 16, symm = FALSE) %>%
    geoCorrection(scl = FALSE)

  # FIXME: rename ID var
  ID <- c(1:ncell(Raster))[!is.na(values(Raster))]
  B <- xyFromCell(Raster, cell = ID)

  # find connections
  connections <- list()
  # for each pair of cells in B
  for (i in 1:nrow(B)) {
    # create a temporal raster for each row with the distance from cell xy to all other cells
    temp <- accCost2(h16, B[i, ])
    index <- which(temp < Dist)
    connections[[i]] <- cbind(ID[i], index, temp[index])
  }

  #Get everything together as a large data frame
  connections <- do.call('rbind', connections)
  connections <- as.data.frame(connections)
  colnames(connections) <- c('from', 'to', 'dist')

  Cons <- list()

  # FIXME: There is a clever tidyverse way to do this...
  for (i in 1:(nlayers(layers_habitat) - 1)) {
    Cons[[i]] <- data.frame(Cell_From = connections$from, Cell_To = connections$to, Time_From = NA, Time_to = NA, Capacity = NA)
    # Why is time zero-indexed here?
    Cons[[i]]$Time_From <- i - 1
    Cons[[i]]$Time_to <- i
    for (j in 1:nrow(connections)) {

      if (
        # connection from is suitable
        (df_habitat %>%
         filter(ID == connections$from[j] & Time == i - 1) %>%
         pull(df_habitat) %>%
         {. == 1})

        &

        # connection to is suitable
        (df_habitat %>%
         filter(ID == connections$to[j] & Time == i) %>%
         pull(df_habitat) %>%
         {. == 1})) {

        # set capacity to one
        Cons[[i]]$Capacity[j] <- 1
      }
    }
    # print layer numbers
    # FIXME: should be informative message
    print(i)
  }

  # convert into one tidy dataframe of only valid edges
  Cons <- Cons %>%
    bind_rows() %>%
    dplyr::filter(Capacity == 1)
  # get cost value of each edge
  Cost <- tibble(
    ID = unique(c(Cons$Cell_To, Cons$Cell_From)),
    cost = layer_cost[ID])

  # add cost information to connection dataframe
  Cons$Cost <- NA
  for (i in 1:nrow(Cons)){
    temp <- Cost %>% dplyr::filter(ID %in% unique(c(Cons$Cell_From[i], Cons$Cell_To[i])))
    if(nrow(temp) == 1){
      Cons$Cost[i] <- temp$cost*2
    }
    if(nrow(temp) == 2){
      Cons$Cost[i] <- temp %>% summarise(cost = sum(cost)) %>% pull(cost)
    }
  }

  ## generate the source and the sink ==========================================

  ForSource <- Cons %>%
    dplyr::filter(Time_From == min(Time_From)) %>%
    mutate(
      Cell_To = Cell_From,
      Cell_From = 0,
      Time_From = -1,
      Time_to = 0,
      Capacity = 1,
      Cost = 0)
  ForSink <- Cons %>%
    dplyr::filter(Time_to == max(Time_to)) %>%
    mutate(
      Cell_From = Cell_To,
      Cell_To = length(values(layers_habitat[[1]][[1]])) + 1,
      Time_From = max(Time_to),
      Time_to = max(Time_to) + 1,
      Capacity = 1,
      Cost = 0)

  Cons <- bind_rows(ForSource, Cons, ForSink) %>% distinct()
  Cons$ID <- 1:nrow(Cons)

  constraintsMatrix <- create_constraints_matrix(Cons, nchains)

  # compute flow solution
  solution <- lpSolve::lp(
    direction = 'min',
    objective.in = Cons[['Cost']],
    const.mat = constraintsMatrix[['lhs']],
    const.dir = constraintsMatrix[['dir']],
    const.rhs = constraintsMatrix[['rhs']])
  # Include solution in edge dataframe
  Cons[['flow']] <- solution[['solution']]

  TestStack <- list()
  for(i in 1:nlayers(layers_habitat)){
    Test <- layer_cost
    values(Test) <- 0
    values(Test)[Cons %>% dplyr::filter(Time_From == (i - 1)) %>% pull(Cell_From) %>%  unique()] <- Cons %>% dplyr::filter(Time_From == (i - 1)) %>% group_by(Cell_From) %>% summarise(flow = sum(flow))  %>% pull(flow)
    TestStack[[i]] <- Test
  }

  TestStack <- do.call('stack', TestStack)
  names(TestStack) <- paste0('T', 1:nlayers(TestStack))

  return(list(Cons = Cons, Stack = TestStack, Soultions = solution))
}

a <- NewFunc(
  layers_habitat = BinSpp[[1]],
  layer_cost = Cost,
  Dist = 400000,
  nchains = 4)

a$Cons %>% group_by(Time_From) %>% summarise(flow = sum(flow))
plot(a$Stack)

b <- NewFunc(Stack = Philips, Dist = 111000, costlayer = PhilCost, nchains = 2)

a$Cons %>%
  group_by(Time_From) %>%
  summarise(flow = sum(flow))
b$Cons %>%
  group_by(Time_From) %>%
  summarise(flow = sum(flow))


NewFuncByNode <- function(Stack = BinSpp[[1]], Dist = 400000, costlayer = Cost, nchains = 8){
  Masklayer <- costlayer
  values(Masklayer) <- ifelse(is.na(values(Masklayer)), NA, 1)
  Stack <- Stack * Masklayer
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

  Suitability <- list()
  for (i in 1:nlayers(Stack)){
    temp <- data.frame(Suitability = values(Stack[[i]]), ID = 1:length(values(Stack[[i]])), Time = i-1)
    Suitability[[i]] <- temp[complete.cases(temp),]
  }
  Suitabilities <- do.call('rbind', Suitability)

  s <- Suitabilities %>% group_by(ID) %>% summarise(SUMA = sum(Suitability)) %>% filter(SUMA > 0)
  Suitabilities <- Suitabilities[Suitabilities$ID %in% s$ID,]

  ###
  Raster <- sum(Stack)

  Raster[values(Raster) > 0] = 1
  Raster[values(Raster) == 0] = NA

  h16 <- transition(Raster, transitionFunction=function(x){1},16,symm=FALSE)

  h16 <- geoCorrection(h16, scl=FALSE)

  ID <- c(1:ncell(Raster))[!is.na(values(Raster))]

  B <- xyFromCell(Raster, cell = ID)

  connections <- list()
  #For each pair of cells in B
  for (i in 1:nrow(B)){
    #Create a temporal raster for each row with the distance from cell xy to all other cells
    temp <- accCost2(h16,B[i,])
    index <- which(temp < Dist)
    connections[[i]] <- cbind(ID[i], index, temp[index])
  }
  #Get everything together as a large data frame
  connections <- do.call('rbind', connections)
  connections <- as.data.frame(connections)
  colnames(connections) <- c('from', 'to', 'dist')

  Cons <- list()

  for(i in 1:(nlayers(Stack)-1)){
    Cons[[i]] <- data.frame(Cell_From = connections$from, Cell_To = connections$to, Time_From = NA, Time_to = NA, Capacity = NA)
    Cons[[i]]$Time_From <- i-1
    Cons[[i]]$Time_to <- i
    for(j in 1:nrow(connections)){
      if((dplyr::filter(Suitabilities, ID == connections$from[j] & Time == i-1) %>% pull(Suitability)) == 1 & (dplyr::filter(Suitabilities, ID == connections$to[j] & Time == i) %>% pull(Suitability)) == 1){
        Cons[[i]]$Capacity[j] <- 1
      }
    }
    print(i)
  }

  Cons <- bind_rows(Cons) %>% dplyr::filter(Capacity == 1)

  Cost <- data.frame(ID = unique(c(c(Cons$Cell_To), c(Cons$Cell_From))), cost = values(costlayer)[unique(c(c(Cons$Cell_To), c(Cons$Cell_From)))])

  Suitabilities <- full_join(Suitabilities, Cost)

  Source_Node <- data.frame(Suitability = nchains, ID = (min(Suitabilities$ID) - 1), Time = (min(Suitabilities$Time) - 1), cost = 0)
  Target_Node <- data.frame(Suitability = nchains, ID = (max(Suitabilities$ID) + 1), Time = (max(Suitabilities$Time) + 1), cost = 0)

  Suitabilities <- bind_rows(Source_Node, Suitabilities, Target_Node) %>% dplyr::filter(Suitability > 0) %>% distinct()

  Suitabilities$TimeSpaceID <- 1:nrow(Suitabilities)

  Cons$Cost <- NA

  for(i in 1:nrow(Cons)){
    temp <- Cost %>% dplyr::filter(ID %in% unique(c(Cons$Cell_From[i], Cons$Cell_To[i])))
    if(nrow(temp) == 1){
      Cons$Cost[i] <- temp$cost*2
    }
    if(nrow(temp) == 2){
      Cons$Cost[i] <- temp %>% summarise(cost = sum(cost)) %>% pull(cost)
    }
  }

  ##Generate the source and the sink

  #Source

  ForSource <- Cons %>% dplyr::filter(Time_From == min(Time_From)) %>% mutate(Cell_To = Cell_From, Cell_From = 0, Time_From = -1, Time_to = 0, Capacity = 1, Cost = 0)

  ForSink <- Cons %>% dplyr::filter(Time_to == max(Time_to)) %>% mutate(Cell_From = Cell_To, Cell_To = length(values(Stack[[1]][[1]])) + 1, Time_From = max(Time_to), Time_to = max(Time_to) + 1, Capacity = 1, Cost = 0)

  Cons <- bind_rows(ForSource, Cons, ForSink) %>% distinct()

  Cons$ID <- 1:nrow(Cons)

  createConstraintsMatrix <- function(edges, total_flow, Nodes) {

    # Edge IDs to be used as names
    names_edges <- edges$ID
    # Number of edges
    numberof_edges <- length(names_edges)

    # Node IDs to be used as names
    names_nodes <- unique(Nodes$TimeSpaceID)
    # Number of nodes
    numberof_nodes <- length(names_nodes)

    # Times id

    name_times <- Nodes$Time %>% unique
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

    # Flow through individual nodes
    constraints$lhs <- Nodes$TimeSpaceID %>%
      length %>%
      diag %>%
      magrittr::set_colnames(Nodes$TimeSpaceID) %>%
      magrittr::set_rownames(Nodes$TimeSpaceID)
    # should be smaller than or equal to
    constraints$dir <- rep('<=', times = nrow(Nodes))
    constraints$dir[1] <- '=='
    constraints$dir[nrow(Nodes)] <- '=='
    # than capacity
    constraints$rhs <- Nodes$Suitability


    #' Build edge flow constraints -----------------------------------------------
    #' For each node, find all edges that go to that node
    #' and all edges that go from that node. The sum of all inputs
    #' and all outputs should be zero. So we set inbound edge coefficients as 1
    #' and outbound coefficients as -1. In any viable solution the result should
    #' be equal to zero.
    #'
    edges$Time_Space_From <-NA
    edges$Time_Space_To <- NA

    for(i in 1:nrow(edges)){
      edges$Time_Space_From[i] <- Nodes %>% dplyr::filter(ID == edges$Cell_From[i], Time == edges$Time_From[i]) %>% pull(TimeSpaceID)
      edges$Time_Space_To[i] <- Nodes %>% dplyr::filter(ID == edges$Cell_To[i], Time == edges$Time_to[i]) %>% pull(TimeSpaceID)
    }

    nodeflow <- matrix(0,
                       ncol = numberof_nodes,
                       nrow = numberof_edges,
                       dimnames = list(names_edges, names_nodes))

    for (i in names_nodes) {
      # input arcs
      edges_in <- edges %>%
        dplyr::filter(Time_Space_To == i) %>%
        dplyr::select(ID) %>%
        unlist
      # output arcs
      edges_out <- edges %>%
        dplyr::filter(Time_Space_From == i) %>%
        dplyr::select(ID) %>%
        unlist

      # output arcs

      # set input coefficients to 1
      nodeflow[rownames(nodeflow) %in% edges_in,
               colnames(nodeflow) == i] <- 1

      # set output coefficients to -1
      nodeflow[rownames(nodeflow) %in% edges_out,
               colnames(nodeflow) == i] <- -1
    }

    # But exclude source and target nodes
    # as the zero-sum flow constraint does not apply to these!
    # Source node is assumed to be the one with the minimum ID number
    # Sink node is assumed to be the one with the maximum ID number
    sourcenode_id <- min(Nodes$TimeSpaceID)
    targetnode_id <- max(Nodes$TimeSpaceID)

    # Exclude them from node flow here
    nodeflow <- nodeflow[nodeflow[,sourcenode_id] ==  0,]
    nodeflow <- nodeflow[nodeflow[,targetnode_id] ==  0,]

    # Add nodeflow to the constraints list
    constraints$lhs <- rbind(constraints$lhs, nodeflow)
    constraints$dir <- c(constraints$dir, rep('==', times = nrow(nodeflow)))
    constraints$rhs <- c(constraints$rhs, rep(0, times = nrow(nodeflow)))


    #' Build initialisation constraints ------------------------------------------
    #' For the source and the target node, we want all outbound nodes and
    #' all inbound nodes to be equal to the sum of flow through the network
    #' respectively
    # Keep node flow values for separate step below
    nodeflow_source <- matrix(nrow = 1, ncol = numberof_nodes, 0,
                              dimnames = list(sourcenode_id, names_nodes))
    nodeflow_source[,sourcenode_id] <- -1


    nodeflow_target <- matrix(nrow = 1, ncol = numberof_nodes, 0,
                              dimnames = list(targetnode_id, names_nodes))
    nodeflow_target[,targetnode_id] <- 1

    # Add initialisation to the constraints list
    constraints$lhs <- rbind(constraints$lhs,
                             source = nodeflow_source,
                             target = nodeflow_target)
    constraints$dir <- c(constraints$dir, rep('==', times = 2))
    # Flow should be negative for source, and positive for target
    constraints$rhs <- c(constraints$rhs, total_flow *-1, total_flow)

    ##Generate the time constraints, find all the nodes en each time

    Timeflow <- matrix(0,
                       nrow = numberof_times,
                       ncol = numberof_nodes,
                       dimnames = list(name_times, names_nodes))

    for (i in name_times) {
      # input arcs
      edges_in <- Nodes %>%
        dplyr::filter(Time == i + 1) %>%
        dplyr::pull(TimeSpaceID)
      # output arcs
      edges_out <- Nodes %>%
        dplyr::filter(Time == i) %>%
        dplyr::pull(TimeSpaceID)

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
    Timeflow <- Timeflow[rownames(Timeflow) != sourcetime_id,]
    Timeflow <- Timeflow[rownames(Timeflow) != targettime_id,]

    constraints$lhs <- rbind(constraints$lhs, Timeflow)
    constraints$dir <- c(constraints$dir, rep('==', times = nrow(Timeflow)))
    constraints$rhs <- c(constraints$rhs, rep(0, times = nrow(Timeflow)))

    ###################
    return(constraints)
  }

  constraintsMatrix <- createConstraintsMatrix(Cons, Nodes = Suitabilities, nchains)

  solution <- lp(
    direction = 'min',
    objective.in = (Cons$Cost),
    const.mat = constraintsMatrix$lhs,
    const.dir = constraintsMatrix$dir,
    const.rhs = constraintsMatrix$rhs)

  # Include solution in edge dataframe
  Cons$flow <- solution$solution

  TestStack <- list()

  for(i in 1:nlayers(Stack)){
    Test<- costlayer
    values(Test) <- 0
    values(Test)[Cons %>% dplyr::filter(Time_From == (i -1)) %>% pull(Cell_From) %>%  unique()] <- Cons %>% dplyr::filter(Time_From == (i - 1)) %>% group_by(Cell_From) %>% summarise(flow = sum(flow))  %>% pull(flow)
    TestStack[[i]] <- Test
  }

  TestStack <- do.call('stack', TestStack)

  names(TestStack) <- paste0('T', 1:nlayers(TestStack))


  return(list(Cons = Cons, Stack = TestStack, Soultions = solution))
}

