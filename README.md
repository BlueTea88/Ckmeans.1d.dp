# Ckmeans.1d.dp.mod

Modified version of the (Ckmeans.1d.dp)[https://github.com/cran/Ckmeans.1d.dp] R package.

Given a vector of `x` values, associated vectors of weights and `z` values for each `x`, and the desired number of clusters, this package aims to cluster the values of `x` such that:

* The weighted sum of squares calculated by the weights and `z` values is minimised; and
* The continuity of `x` values is preserved (only neighbouring values can be grouped together)
