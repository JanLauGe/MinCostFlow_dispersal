#' Edges for the minimum cost flow problem of Steven J Phillips
#'
#' A dataset containing the graph problem outlined in "Target-based site
#' prioritization under climate change" by Steven J Phillips (2017). The graph
#' is defined by a data frame of edge information.
#'
#' @format A data frame with 14 rows and 5 columns:
#' \describe{
#'   \item{edge_id}{primary key, ID of each edge}
#'   \item{node_from}{foreign key, ID of the node that the edge originates at}
#'   \item{node_to}{foreign key, ID of the node that the edge connects to}
#'   \item{edge_capacity}{capacity for flow through the given edge}
#'   \item{edge_cost}{cost of allocating flow through the given edge}
#'   ...
#' }
#' @source \url{https://github.com/JanLauGe/MinCostFlow_dispersal/blob/master/literature/Phillips_2017_targetbased_site_prioritization_under_climate_change.pdf}
"phillips_problem"
