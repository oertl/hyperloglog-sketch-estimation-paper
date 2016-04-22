# Maximum likelihood estimation of cardinalities from HyperLogLog sketches
[Paper](https://github.com/oertl/hyperloglog-sketch-estimation-paper/blob/master/paper/paper.pdf) about the estimation of cardinalities from HyperLogLog sketches. (Work is still in progress.)

## Abstract
This paper presents a new method to estimate the cardinality of a multiset recorded by HyperLogLog sketches. The maximum likelihood method is applied to the probability mass function under the Poisson model. A numerical algorithm is presented that is able to calculate the estimate efficiently. Simulations finally show that this new estimation procedure works well and is inherently unbiased over a wide range of cardinalities.

## Key results
The following two figures show how the new algorithm improves cardinality estimation for a HyperLogLog sketch with 4096 registers using 32bit hash values as input.
* The relative estimation error for the cardinality when using the original HyperLogLog algorithm:
![original](/paper/original_estimate.png)
* The relative estimation error when using the new estimation algorithm based on maximumm likelihood estimation under the Poisson  model:
![original](/paper/max_likelihood_estimate_12_20.png)


