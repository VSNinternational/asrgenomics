geno.apple2 <- geno.apple[1:20, 1:50]
test_that("G calculation works", {
  G <- G.matrix(M = geno.apple2, method = "VanRaden", na.string="NA")$G
  expect_equal(
    G[1:2,1:2],
    matrix(c(1.14002041, -0.09223346, -0.09223346, 1.05094182),
           nrow = 2,
           dimnames = list(c("A325-1", "A325-2"), c("A325-1", "A325-2"))
    )
  )
  Gs <- G.matrix(M = geno.apple2, method = "VanRaden", na.string="NA", sparseform = TRUE, message = FALSE)$G
  expect_equal(
    sparse2full(Gs),
    {attributes(G)$rowNames <- rownames(G) -> attributes(G)$colNames ; G}
  )
})
test_that("traps work", {
  geno.applewr <- as.data.frame(geno.apple)
  expect_error(G.matrix(M = geno.applewr, method = "VanRaden", na.string="NA"))
  geno.applewr <- geno.apple2
  rownames(geno.applewr) <- c()
  expect_error(G.matrix(M = geno.applewr, method = "VanRaden", na.string="NA"))
  geno.applewr <- geno.apple2
  colnames(geno.applewr) <- c()
  expect_error(G.matrix(M = geno.applewr, method = "VanRaden", na.string="NA"))
  geno.applewr <- geno.apple2[, 1:10]
  expect_warning(G.matrix(M = geno.applewr, method = "VanRaden", na.string="NA"))
})
