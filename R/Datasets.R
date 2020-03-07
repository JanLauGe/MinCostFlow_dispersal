
#' A List with two Stacks each for an hypothetic species
#'
#' A dataset containing two stacks in a list, each contains eight rasters, each represents
#' a timeslice of the distribution of the species
#'
#' @format A list with two stacks of the distribution of two ypothetical species:
#' \describe{
#'   \item{SeciesA}{A stack with eight time slices for species A}
#'   \item{SpeciesB}{A stack with eight time slices for species b}
#' }
#'
"BinSpp"

#' A raster with a simulated cost for our hypothetical species
#'
#' This raster has areas with cost Zero (protected areas), unavailable cells
#' (cities or agricultural areas) and cells with cost.
#'
#' @format Raster
#' \describe{
#'   \item{Cost}{simulated cost layer}
#'}
"Cost"

#' A rasterstack with two timeslices of an hypothetical species
#'
#' A dataset containing a stack with two rasters, each represents
#' a timeslice of the distribution of the species
#'
#' @format A rasterstacks of the distribution an hypothetical species:
#' \describe{
#'   \item{Phillips}{A stack with two time slices for a Species}
#' }
#'
"Phillips"

#' A raster with a simulated cost for an hypothetical species
#'
#' This raster has areas with cost Zero (protected areas), unavailable cells
#' (cities or agricultural areas) and cells with cost.
#'
#' @format Raster
#' \describe{
#'   \item{PhilCost}{simulated cost layer}
#'}
"PhilCost"
