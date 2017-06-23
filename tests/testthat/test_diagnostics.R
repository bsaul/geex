context("Geex diagnostics")

test_eefun1 <- function(data){
  function(theta){
    with(data,
         c(Y1 - theta[1],
           (Y1 - theta[1])^2 - theta[2] )
    )
  }
}

test_that("check_GFUN function works", {
  myList <- list(eeFUN = test_eefun1, splitdt = split(geexex, 1:nrow(geexex)))

  temp <- estimate_equations(
    test_eefun1,
    data = geexex,
    roots = c(2,2)
  )

  test <- check_GFUN(myList, temp$parameters)

  expect_equal(length(test), 2)

  expect_equal(test, c(0, 0), tolerance = 1e-5, check.attributes = FALSE)

})



