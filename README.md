# New cardinality estimation algorithms for HyperLogLog sketches
[Paper](http://oertl.github.io/hyperloglog-sketch-estimation-paper/paper.pdf) about the estimation of cardinalities from HyperLogLog sketches. (Work is still in progress.)

## Abstract
This paper presents new methods to estimate the cardinalities of multisets recorded by HyperLogLog sketches. A theoretically founded extension to the original estimate is presented that eliminates the bias for small and large cardinalities. Based on the maximum likelihood principle another unbiased method is derived together with a robust and efficient numerical algorithm to calculate the estimate. The maximum-likelihood method is also appropriate to improve cardinality estimates of set intersections compared to the inclusion-exclusion principle. The new methods are demonstrated and verified by extensive simulations.

## Key results
The following two figures show how the new algorithm improves cardinality estimation for a HyperLogLog sketch with 4096 registers using 32-bit hash values as input.
* The relative estimation error for the cardinality when using the original HyperLogLog algorithm:
![original algorithm relative error chart](https://github.com/oertl/hyperloglog-sketch-estimation-paper/raw/master/paper/original_estimate.png)
* The relative estimation error when using the new estimation algorithm based on maximumm likelihood estimation under the Poisson  model:
![new algorithm relative error chart](https://github.com/oertl/hyperloglog-sketch-estimation-paper/raw/master/paper/max_likelihood_estimate_12_20.png)

The new algorithm is also fast. The following figure shows the average computation time for the cardinality estimate from HyperLogLog sketches with 4096 registers working with a 32-bit (q=20) or 64-bit (q=52) hash function, respectively.
![performance chart](https://github.com/oertl/hyperloglog-sketch-estimation-paper/raw/master/paper/max_likelihood_avg_exec_time.png)
