test_that("selectX() from SelectSim genrates same results on LUAD dataset", {

# Load the dataset
data(luad_run_data, package = "SelectSim")
data(luad_result, package = "SelectSim") 
# result of applying function
  outdat <- selectX(  M = luad_run_data$M,
                        sample.class = luad_run_data$sample.class,
                        alteration.class = luad_run_data$alteration.class,
                        n.cores = 1,
                        min.freq = 10,
                        n.permut = 1000,
                        lambda = 0.3,
                        tao = 1,
                        save.object = FALSE,
                        verbose = FALSE,
                        estimate_pairwise = FALSE,
                        maxFDR = 0.25)

# Check if the result dataframe is same
  expect_true(identical(dim(luad_result), dim(outdat$result)))
# Check if the number of signficant pairs is the same
  true_p <- nrow(luad_result[luad_result$FDR==TRUE,])
  test_p <- nrow(outdat$result[outdat$result$FDR==TRUE,])
  expect_equal(true_p,test_p)
})
#> Test passed 😀