## ---- echo = FALSE, message = FALSE, warning=FALSE-----------------------
library(geex)

## ----SB6_eefun, echo = TRUE----------------------------------------------
eefun <- function(data, k = 1.5){
  function(theta){
    x <- data$Y1 - theta[1]
    if(abs(x) <= k) x else sign(x) * k
  }
}

## ----multiroot_example, echo = TRUE, message=FALSE-----------------------
multiroot_results <- estimate_equations(
  eeFUN = eefun, 
  data  = geexex,
  rootFUN_control = list(start = 3))

## ----uniroot_example, echo = TRUE, message=FALSE-------------------------
uniroot_results <- estimate_equations(
  eeFUN = eefun, 
  data  = geexex,
  rootFUN = stats::uniroot,
  rootFUN_control = list(interval = c(0, 10)))

## ----compare_results-----------------------------------------------------
multiroot_results$estimates - uniroot_results$estimates

