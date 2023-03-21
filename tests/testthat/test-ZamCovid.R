context("ZamCovid")


test_that("can run ZamCovid model", {

  p <- ZamCovid_parameters(start_date = "2020-02-01")
  mod <- ZamCovid$new(p, 0, 5, seed = 1L)
  info <- mod$info()

  y <- mod$run(400)

  expect_equal(c(120, 5), dim(y))
})


test_that("can run basic model", {

  p <- basic_parameters(model_end = 200)
  mod <- basic$new(p, 0, 5)
  info <- mod$info()

  y <- mod$run(200)

  expect_equal(c(104, 5), dim(y))
})
