# dependencies
if (!require('pacman')) install.packages('pacman')
if (!require('devtools')) install.packages('devtools')
if (!require('QuadCostAmpl')) devtools::install_github('derek-corcoran-barrios/QuadCostAmpl')
pacman::p_load(devtools, tidyverse, magrittr, gdistance, raster, lpSolve, QuadCostAmpl, tictoc)

source('R/funs_min_cost_flow.R')

# example data
data('BinSpp')
data('Cost')

# let's have a look
layers_habitat = BinSpp[[1]]
#layers_habitat <- get_dummy_data()[[1]]
layer_cost = Cost
#layer_cost <- get_dummy_data()[[2]]
Dist = 240000
nchains = 4




connections <- edge_distance_limit(layer_cost = layer_cost, layers_habitat = layers_habitat, Dist = Dist)



#########################
#########################

# data prep --------------------------------------------------------------------
# filter habitat suitability layer for sites with dispersal cost
layers_habitat <- layers_habitat * !is.na(layer_cost)

# create a long list of cells across time slices
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

### Generate a DF that keeps the ID of the raster cell and the new cell ID

df_IDs <- df_habitat %>% rename(node_to = cell_id) %>% dplyr::select(node_to, raster_id)

# extract cost value of each edge
# TODO: there must be a better way to do this without replicating data * n_timeslices
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

# get all valid edge connections,
# i.e. edges between two suitable cells in two adjacent time slices
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
# add edge id and capacity


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


# add source and target edges
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

# sovler -----------------------------------------------------------------------
# create constraints matrix
constraints_matrix <- create_constraints_matrix(
  edges = edges_formatted,
  total_flow = nchains)

# run lpSolve to find best solution
solution <- lp(
  direction = 'min',
  objective.in = edges_formatted[['edge_cost']],
  const.mat = constraints_matrix[['lhs']],
  const.dir = constraints_matrix[['dir']],
  const.rhs = constraints_matrix[['rhs']])
solution


# visualize --------------------------------------------------------------------
# combine with edge information
edges_solved <- bind_cols(edges_formatted, solution = solution[['solution']]) %>% left_join(df_IDs) %>% dplyr::filter(solution != 0)


TestStack <- list()


for(i in 1:nlayers(layers_habitat)){
  Test<- layer_cost
  values(Test) <- 0
  values(Test)[edges_solved %>% dplyr::filter(timeslice_to == (i)) %>% pull(raster_id) %>%  unique()] <- edges_solved %>% dplyr::filter(timeslice_to == (i)) %>% group_by(raster_id) %>% summarise(flow = sum(solution))  %>% pull(flow)
  TestStack[[i]] <- Test
}

TestStack <- do.call("stack", TestStack)

names(TestStack) <- paste0("T", 1:nlayers(TestStack))

plot(TestStack)
plot(max(TestStack))
#####
#########

#st_as_sf(as(extent(Cost), 'SpatialPolygons'))

