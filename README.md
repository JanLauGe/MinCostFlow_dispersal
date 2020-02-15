# MinCostFlow algorithm for species dispersal

## Project purpose
"The idea is to make Minimum Cost Network Flow work through time, so there are more layers since the flow conservation has to apply not only through the network and space but also through time."

## Inputs
1. rasters of habitat suitability at each time slice `t`.
1. raster of cost of dispersal to(?) each site.
1. Dist(?)
1. nchains(number of dispersal events?)

## Algorithm
Find the relative number pf  of dispersal flow paths going through a particular site?

## Outputs
A list of sites and their importance for network dispersal, ranging from `0` to `1`.

## References
[Research note on MinCostFlow for dispersal chains](literature/Phillips_2017_targetbased_site_prioritization_under_climate_change.pdf)
