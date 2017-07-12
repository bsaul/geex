# Overview

`geex` provides an extensible API for estimating parameters and their covariance from a set of estimating functions (M-estimation). M-estimation theory has a long history (see the [M-estimation bibliography](https://bsaul/github.io/geex/docs/articles/articles/mestimation_bib.html)). For an excellent introduction, see the primer by L.A. Stefanski and D.D. Boos,  "The Calculus of M-estimation" ([The American Statistician (2002), 56(1), 29-38)](http://www.jstor.org/stable/3087324?seq=1#page_scan_tab_contents); [also available here](http://www4.stat.ncsu.edu/~boos/papers/mest6.pdf)). 

M-estimation encompasses a broad swath of statistical estimators and ideas including:

* the empirical "sandwich" variance estimator
* generalized estimating equations (GEE)
* many maximum likelihood estimators
* robust regression
* and many more

`geex` can implement all of these using a user-defined estimating function. 

> If you can specify a set of unbiased estimating equations, `geex` can do the rest.

# Installation

To install the current version:

```
devtools::install_github("bsaul/geex")
```

# Usage

Start with the examples in the [package introduction](https://bsaul.github.io/geex/articles/00_geex_intro.html) (also accessible in R by `vignette('00_geex_intro')`). 
