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
  mutate(cell_id = cell_id + max(cell_id) * time)

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
    timeslice_to = .x)) %>%
  # add edge id and capacity
  mutate(
    edge_id = row_number(),
    edge_capacity = 1) %>%
  # add cost
  left_join(df_cost, by = c('node_to' = 'cell_id'))

#' - edge_id (primary key)
#' - node_from (foreign key)
#' - node_to (foreign key)
#' - edge_capacity (int)
#' - edge_cost (int)
