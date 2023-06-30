context("Geex diagnostics")

test_eefun1 <- function(data) {
  function(theta) {
    with(
      data,
      c(
        Y1 - theta[1],
        (Y1 - theta[1])^2 - theta[2]
      )
    )
  }
}

test_that("check_GFUN function works", {
  temp <- m_estimate(
    estFUN = test_eefun1,
    data = geexex,
    root_control = new("root_control", .options = list(start = c(2, 2)))
  )

  expect_silent({
    gtest <- grab_GFUN(temp@basis)
  })
  expect_is(gtest, "function")

  test <- diagnose_roots(gtest, temp@estimates)
  expect_equal(length(test), 2)
  expect_equal(test, c(0, 0), tolerance = 1e-5, check.attributes = FALSE)
})
