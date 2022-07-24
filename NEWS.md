# geex 1.1.1

* Removes `ivpack` dependency from additional examples vignette
in order to get geex back on CRAN.
* Removes the internal `compute_sum_of_list` and `compute_pairwise_sum_of_list`
functions from exports.

# geex 1.1.0

* adds a `Asolver` argument to `m_estimate` to allow for using a different to compute the inverse of the A ("bread") matrix.

# geex 1.0.11

* renames the `geex-summary` class to `geex_summary`
* adds `ee_i` (observed estimating function values) to the component object
* various other bug fixes

# geex 1.0.5

* Added `grab_design_levels` to find all factor-levels in a design matrix. This output can then be passed to the `xlev` arg in `model.matrix` or more importantly through `grab_psiFun` for the models. This fixes issue 41. 
   - Notably, if this seems to work well, then perhaps `grab_design_levels` could be added to the body of `grab_psiFun` to automatically handle this issue.
* Edited `binomial_integrand` to hopefully handle singleton groups better. It will now attempt to coerce back to a one-row matrix before calling `apply`.

# geex 1.0.4

* fixes issue when using `grab_psiFUN` for `glm` objects and non-grouped data.
* removes option to use weights with `grab_psiFUN`, as it was not clear that was implemented correctly.

# geex 1.0.3

* adds basic examples for most functions, though to be sure, the vignettes provide more useful examples.

# geex 1.0.2

* requires R >= 3.3

# geex 1.0.1

* add a `call` slot to the S4 `geex` object. Now the `update` function can be used to update elements of an `m_estimate` call.

# geex 1.0.0

* implements an `S4` system throughout `geex`
* `estimate_equations` becomes `m_estimate`. See documentation for changes to arguments. Notably, `eeFUN` becomes `estFUN`
* `make_eeFUN` functions become `grab_psiFUN`
* speeds up the summation of list of matrices with `compute_sum_of_list`
* plus many more updates and vignettes

# geex 0.3.0

* adds a `weights` argument to `estimate_equations` for faster computations with grouped data. See the [weights vignette](https://bsaul.github.io/geex/articles/v04_weights.html) for a demonstration.

# geex 0.2.2

* changes names of list items used in a `corrections_list` in `estimate_equations`. Each item of the `correction_list` must itself be a list at least one item: `correctFUN`. Additional arguments to `correctFUN` may be passed via `correctFUN_control`.

# geex 0.2.1

* adds a vignette explaining how to use a different root finding algorithm via `rootFUN` argument in `estimate_equations`
* fixes issues where different `rootFUN`s would not work:
    * `roots` argument no longer needs to be set if `compute_roots = TRUE`. Instead, starting values are passed via the `rootFUN_control` list. 
    * Removes the `start` argument from `compute_eeroots`; set this option in `rootFUN_control`.

# geex 0.2.0

* overhauls which arguments are passed to `estimate_equations` and how these arguments are parsed internally. See this function's documentation for details
* adds `geexex` dataset for use in examples

# geex 0.1.0

* the inital `geex` release
