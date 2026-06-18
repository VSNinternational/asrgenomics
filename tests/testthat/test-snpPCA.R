test_that("SNP PCA works", {
  SNP_pca <- snp.pca(M = geno.apple, ncp = 10)
  expect_equal(ncol(SNP_pca$pca.scores), 10)
  expect_equal(nrow(SNP_pca$eigenvalues), 10)
  grp <- as.factor(pheno.apple$Family)
  expect_no_error(
    SNP_pca_grp <- snp.pca(M = geno.apple, groups = grp, label = TRUE, ellipses = TRUE)
  )
  expect_no_error(
    suppressWarnings(SNP_pca_grp <- snp.pca(M = geno.apple, groups = grp, label = FALSE))
  )
  expect_no_error(
    suppressWarnings(SNP_pca_grp <- snp.pca(M = geno.apple, label = TRUE))
  )
})
test_that("traps work", {
  expect_error(
    snp.pca(M = as.data.frame(geno.apple), ncp = 10)
  )
  geno.applewr <- geno.apple
  rownames(geno.applewr) <- c()
  expect_error(
    snp.pca(M = geno.applewr, ncp = 10)
  )
  geno.applewr <- geno.apple
  colnames(geno.applewr) <- c()
  expect_error(
    snp.pca(M = geno.applewr, ncp = 10)
  )
  geno.applewr <- geno.apple
  geno.applewr[1, 1] <- NA
  expect_error(
    snp.pca(M = geno.applewr, ncp = 10)
  )
  expect_error(
    snp.pca(M = geno.apple, ncp = 2829)
  )
})
