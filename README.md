# sfclust: Bayesian Spatial Functional Clustering

**sfclust** provides a Bayesian framework for clustering spatio-temporal data,
supporting both Gaussian and non-Gaussian responses. The approach enforces spatial
adjacency constraints, ensuring that clusters consist of neighboring regions with
similar temporal dynamics.

The package implements the methodology described in *"Bayesian Spatial Functional
Data Clustering: Applications in Disease Surveillance"*, available at
[arXiv:2407.12633](https://arxiv.org/abs/2407.12633). In addition to the core
clustering algorithm, `sfclust` offers tools for model diagnostics, visualization,
and result summarization.

## Installation

`sfclust` relies on the [`INLA`](https://www.r-inla.org/download-install) package for
efficient Bayesian inference. Install it with:

```r
install.packages("INLA", dependencies = TRUE,
  repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable")
)
```

Once `INLA` is installed, you can install the development version of `sfclust` from
GitHub:

```r
devtools::install_github("ErickChacon/sfclust")
```

## Basic Usage

Suppose you have a spatio-temporal `stars` object named `stars_object` that contains
variables such as `cases` and `expected` (the expected number of cases). The
following code fits a spatial functional clustering model, where each cluster’s mean
trend is modeled with a temporal random walk and an unstructured random effect:

```r
form <- cases ~ f(idt, model = "rw1") + f(id, model = "iid")
result <- sfclust(stars_object, formula = form, family = "poisson", E = expected, niter = 1000)
result
```

## Acknowledgments

We thank the authors of *"Bayesian Clustering of Spatial Functional Data with
Application to a Human Mobility Study During COVID-19"*, by Bohai Zhang, Huiyan Sang,
Zhao Tang Luo, and Hui Huang (*Annals of Applied Statistics*, 2023), for making their
supplementary code publicly available
([DOI:10.1214/22-AOAS1643SUPPB](https://doi.org/10.1214/22-AOAS1643SUPPB)). Our
implementation builds upon their clustering algorithm and uses their code for
generating spanning trees. We are grateful for their contributions and inspiration.
