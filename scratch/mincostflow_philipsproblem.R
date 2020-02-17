# dependencies
if (!require('pacman')) install.packages('pacman')
if (!require('devtools')) install.packages('devtools')
if (!require('QuadCostAmpl')) devtools::install_github('derek-corcoran-barrios/QuadCostAmpl')
pacman::p_load(devtools, tidyverse, magrittr, gdistance, raster, lpSolve, QuadCostAmpl, igraph, tictoc)
source('functions.R')

# example data -----------------------------------------------------------------
# get edge data for phillips problem
edges_phillips <- data_phillips_problem()
# create constraints matrix
constraints_matrix <- create_constraints_matrix(
  edges = edges_phillips,
  total_flow = 2)

# sovler -----------------------------------------------------------------------
# run lpSolve to find best solution
solution <- lp(
  direction = 'min',
  objective.in = edges_phillips[['edge_cost']],
  const.mat = constraints_matrix[['lhs']],
  const.dir = constraints_matrix[['dir']],
  const.rhs = constraints_matrix[['rhs']])
# print vector of flow by edge
solution[['solution']]

# visualize output -------------------------------------------------------------
# make edgelist into graph object
g <- graph_from_edgelist(as.matrix(edges_phillips[,c('node_from','node_to')]))
# add properties
E(g)$capacity <- edges_phillips$edge_capacity
E(g)$cost <- edges_phillips$edge_cost
# get some colours in to visualise cost
E(g)$color[E(g)$cost == 0] <- 'royalblue'
E(g)$color[E(g)$cost >= 1] <- 'firebrick'

# visualize graph
plot(g, edge.label = E(g)$capacity)
plot(g, edge.label = E(g)$cost)

# visualize solution:
# flow as edge size,
# cost as colour
E(g)$flow <- solution$solution
plot(g, edge.width = E(g)$flow * 10)


