# What does `geex` do?

`geex` provides a framework for estimating parameters from a set of estimating equations (M-estimation). M-estimation theory goes back to Huber and others, but L.A. Stefanski and D.D. Boos provide an excellent introduction in "The Calculus of M-estimation" ([The American Statistician (2002), 56(1), 29-38)](http://www.jstor.org/stable/3087324?seq=1#page_scan_tab_contents)). Once a set of estimating equations, such the score equations in maximum likelihood, are definted, M-estimation provides the theoretical crank to turn that yields asymptotically Normal estimates. 

`geex` takes a user-defined set of estimating equations and returns parameter and covariance estimates. The covariance estimates are the usual "sandwich" estimates. These are known to underestimate the true variance in small samples, so the package also computes finite sample corrections proposed by M.P. Fay, and B.I. Graubard in "Small-Sample adjustments for wald-type tests using sandwich estimators" (Biometrics (2001), 57(4), 1198-1206).

# How does `geex` work?

The key piece is the user-defined set of estimating equations, which is a function that takes in group-level data and returns a function with respect to the parameters.

For example, to estimate the mean and variance, one can solve the following set of estimating equations:

```
example_ee <- function(data){
    function(theta){
      with(data,
        c(Y - theta[1],
         (Y - theta[1])^2 - theta[2] )) 
    }
}
```

`example_ee` is then passed to `geex::estimate_equations()` for computations.

The best way to learn `geex` is to look at the examples in the vignette:

```
vignette('geex_intro')
```

This file walks through all the examples in Stefanski and Boos. 
