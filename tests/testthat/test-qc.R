geno <- matrix(
  c(
    -9, -9, 1, 1, -9, -9,
    -9, -9, 0, 1, 1, 0,
    -9, 1, 2, 1, 1, 1,
    -9, 1, 2, 1, 1, 1,
    -9, 2, 2, 1, 2, 2
  ),
  nrow = 5,
  byrow = TRUE,
  dimnames = list(
    c("I1", "I2", "I3", "I4", "I5"),
    c("M1", "M2", "M3", "M4", "M5", "M6")
  )
)
dummymap <- dummy.map_(marker.id = colnames(geno), message = FALSE)
test_that("filters work", {
  silent_(
    M_filter <-
      qc.filtering(
        M = geno,
        map = dummymap,
        chrom = "chrom", marker = "marker", pos = "pos",
        marker.callrate = 0.39,
        ind.callrate = 0.49,
        maf = 0.26,
        heterozygosity = .99999,
        Fis = .45,
        impute = TRUE,
        na.string = -9,
        message = TRUE
      )
  )
  expect_equal(
    M_filter$M.clean,
    geno[2:5, 6, drop = FALSE]
  )
  expect_equal(
    M_filter$map,
    dummymap[6, ]
  )
})
test_that("imputation works", {
  silent_(
    M_filter <-
      qc.filtering(
        M = geno,
        impute = TRUE,
        na.string = -9,
        message = TRUE
      )
  )
  expect_equal(M_filter$M.clean[1, 1], 1.33)
})
test_that("it understands NA strings", {
  genomod <- geno
  genomod[genomod == -9] <- NA
  expect_no_error(
    M_filter <-
      qc.filtering(
        M = genomod,
        na.string = NA,
        message = FALSE
      )
  )
  expect_no_error(
    M_filter <-
      qc.filtering(
        M = genomod,
        na.string = "NA",
        message = FALSE
      )
  )
})
test_that("it only accepts 0, 1, 2, and NA", {
  genomod <- geno
  genomod[genomod == -9] <- 10
  expect_error(
    M_filter <-
      qc.filtering(
        M = genomod,
        message = FALSE
      )
  )
})
test_that("plots can be ommited", {
  genomod <- geno
  genomod[genomod == -9] <- 10
  M_filter <-
    qc.filtering(
      M = geno,
      na.string = -9,
      plots = FALSE,
      message = FALSE
    )
  expect_null(M_filter$plot.missing.ind)
  expect_null(M_filter$plot.missing.SNP)
  expect_null(M_filter$plot.heteroz)
  expect_null(M_filter$plot.Fis)
  expect_null(M_filter$plot.maf)
})
test_that("traps work", {
  dummymapwr <- dummymap
  dummymapwr$marker[1] <- "null"
  expect_error(
    M_filter <-
      qc.filtering(
        M = geno,
        map = dummymapwr,
        chrom = "chrom", marker = "marker", pos = "pos",
        message = FALSE
      )
  )
  expect_error(
    M_filter <-
      qc.filtering(
        M = geno,
        map = dummymap,
        chrom = "chromo", marker = "marker", pos = "pos",
        message = FALSE
      )
  )
  expect_error(
    M_filter <-
      qc.filtering(
        M = geno,
        map = dummymap,
        message = FALSE
      )
  )
  expect_error(
    M_filter <-
      qc.filtering(
        M = geno,
        map = dummymap,
        heterozygosity = -1,
        message = FALSE
      )
  )
  expect_error(
    M_filter <-
      qc.filtering(
        M = geno,
        map = dummymap,
        marker.callrate = 0.39,
        ind.callrate = 0.49,
        maf = 0.26,
        heterozygosity = .99999,
        Fis = -1,
        impute = TRUE,
        na.string = -9,
        message = FALSE
      )
  )
  expect_error(
    M_filter <-
      qc.filtering(
        M = geno,
        map = dummymap,
        maf = -1,
        message = FALSE
      )
  )
  expect_error(
    M_filter <-
      qc.filtering(
        M = geno,
        map = dummymap,
        ind.callrate = -1,
        message = FALSE
      )
  )
  expect_error(
    M_filter <-
      qc.filtering(
        M = geno,
        map = dummymap,
        marker.callrate = -1,
        message = FALSE
      )
  )
  genowr <- geno
  rownames(genowr) <- NULL
  expect_error(
    M_filter <-
      qc.filtering(
        M = genowr,
        map = dummymap,
        message = FALSE
      )
  )
  genowr <- geno
  colnames(genowr) <- NULL
  expect_error(
    M_filter <-
      qc.filtering(
        M = genowr,
        map = dummymap,
        message = FALSE
      )
  )
  expect_error(
    M_filter <-
      qc.filtering(
        M = as.data.frame(genowr),
        map = dummymap,
        message = FALSE
      )
  )
  expect_error(
    M_filter <-
      qc.filtering(
        M = as.data.frame(genowr),
        map = dummymap,
        ref = c("a"),
        message = FALSE
      )
  )
})
