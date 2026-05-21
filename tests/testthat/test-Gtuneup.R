G <- G.matrix(
  M = geno.apple[1:20, 1:50],
  method = "VanRaden",
  na.string = "NA"
)$G
dummymat <-
  matrix(
    c(1, 0, 0, 0, 1, .25, 0, .25, 1),
    byrow = TRUE,
    nrow = 3,
    dimnames = list(1:3, 1:3)
  )
A <- diag(nrow = nrow(dummymat))
A[lower.tri(A)] <- A[upper.tri(A)] <- runif(3, min = 0, max = .5)
rownames(A) <- colnames(A) <- rownames(dummymat)
test_that("tuneup works", {
  Gb <- G.tuneup(G = dummymat, blend = TRUE, pblend = .1)$Gb
  expect_equal(
    Gb,
    (dummymat * .9) + diag(3) * .1
  )
  Gb <- G.tuneup(G = G, bend = TRUE)
  expect_gt(Gb$rcnb, Gb$rcn0)
  Gb <- G.tuneup(G = dummymat, A = A, align = TRUE)$Gb
  expect_false(identical(Gb, dummymat))
  Gb <- G.tuneup(G = dummymat, A = A, blend = TRUE)$Gb
  expect_false(identical(Gb, dummymat))
  Gb <- G.tuneup(G = dummymat, A = A, blend = TRUE, sparseform = TRUE)$Gb.sparse
  expect_equal(ncol(Gb), 3)
})
test_that("traps work", {
  expect_error(
    G.tuneup(G = dummymat, align = TRUE)$Gb
  )
  expect_warning(
    G.tuneup(G = dummymat, A = A[1:2, 1:2], align = TRUE)$Gb
  )
  expect_error(
    suppressWarnings(
      G.tuneup(G = dummymat, A = A[1:2, 1:2], blend = TRUE)$Gb
    )
  )
  expect_error(
    G.tuneup(G = dummymat, pblend = -1)$Gb
  )
  largemat <- diag(1501, nrow = 1501)
  dimnames(largemat) <- list(1:1501, 1:1501)
  Gb <- G.tuneup(G = largemat, blend = TRUE, determinant = F)
  expect_null(Gb$det0)
  Gb <- G.tuneup(G = largemat, blend = TRUE, determinant = T)
  expect_false(is.null(Gb$det0))
  expect_error(
    G.tuneup(G = dummymat, A = A[1:2, ], align = TRUE)$Gb
  )
  Awr <- A
  colnames(Awr) <- c()
  expect_error(
    G.tuneup(G = dummymat, A = Awr, align = TRUE)$Gb
  )
  Awr <- A
  rownames(Awr) <- c()
  expect_error(
    G.tuneup(G = dummymat, A = Awr, align = TRUE)$Gb
  )
  expect_error(
    G.tuneup(G = dummymat, A = as.data.frame(A), align = TRUE)$Gb
  )
  expect_error(
    G.tuneup(G = dummymat[1:2, ], blend = TRUE)$Gb
  )
  expect_error(
    G.tuneup(G = dummymat)$Gb
  )
  expect_error(
    G.tuneup(G = dummymat, blend = TRUE, bend = TRUE, align = TRUE)$Gb
  )
  expect_error(
    G.tuneup(G = dummymat, blend = TRUE, eig.tol = -1)$Gb
  )
  expect_error(
    G.tuneup(G = dummymat, blend = TRUE, pblend = -1)$Gb
  )
  expect_error(
    G.tuneup(G = dummymat[1:2, ], blend = TRUE)$Gb
  )
  dummymatwr <- dummymat
  colnames(dummymatwr) <- c()
  expect_error(
    G.tuneup(G = dummymatwr, blend = TRUE)$Gb
  )
  dummymatwr <- dummymat
  rownames(dummymatwr) <- c()
  expect_error(
    G.tuneup(G = dummymatwr, blend = TRUE)$Gb
  )
  expect_error(
    G.tuneup(G = as.data.frame(dummymat), blend = TRUE)$Gb
  )
  Awr <- A
  colnames(Awr) <- rownames(Awr) <- 4:6
  expect_error(
    suppressWarnings(G.tuneup(G = dummymat, A = Awr, blend = TRUE)$Gb)
  )
})
