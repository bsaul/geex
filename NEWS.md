# geex 0.3.0

* adds a `weights` argument to `estimate_equations` for faster computations with grouped data. See the [weights vignette](https://bsaul.github.io/geex/docs/articles/weights.html) for a demonstration.

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
