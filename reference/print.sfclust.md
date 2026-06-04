# Print method for sfclust objects

Prints details of an sfclust object, including the (i) within-cluster
formula; (ii) hyperparameters used for the MCMC sample such as the
number of clusters penalty (q) and the movement probabilities
(move_prob); (iii) the number of movement type counts during the MCMC
sampling; and (iv) the log marginal likelihood of the model of the last
clustering sample.

## Usage

``` r
# S3 method for class 'sfclust'
print(x, ...)
```

## Arguments

- x:

  An object of class 'sfclust'.

- ...:

  Additional arguments passed to `print.default`.

## Value

Invisibly returns the input `sfclust` object `x`. The function also
prints a summary of:

- the within-cluster model formula,

- clustering hyperparameters,

- movement counts from the MCMC sampler,

- and the log marginal likelihood of the selected sample.
