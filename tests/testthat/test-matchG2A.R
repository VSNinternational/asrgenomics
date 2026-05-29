dummymat <-
  matrix(c(1,0,0,0,
           0,1,.25,0,
           0,.25,1,0,
           0, 0, 0, 1),
         byrow = TRUE, nrow = 4,
         dimnames = list(1:4, 1:4))
A <- diag(nrow = nrow(dummymat))
A[lower.tri(A)] <- A[upper.tri(A)] <- c(.5,.25,0, 0, 0, 0)
rownames(A) <- colnames(A) <- rownames(dummymat)
test_that("mathing G2A works", {
  check <- match.G2A(A = A, G = dummymat, clean = TRUE, ord = TRUE, mism = TRUE,
                   RMdiff = TRUE)
  expect_equal(check$Aclean, A)
  check <- match.G2A(A = A, G = dummymat, clean = TRUE, ord = TRUE, mism = TRUE,
                   RMdiff = FALSE)
  expect_null(check$RM)
  check <- match.G2A(A = A, G = dummymat[1:2, 1:2], clean = TRUE, ord = FALSE, mism = FALSE,
                   RMdiff = FALSE)
  expect_equal(check$Aclean, A[1:2,1:2])
  expect_message(
    check <- match.G2A(A = A[3:1, 3:1], G = dummymat, ord = FALSE)
  )
  expect_message(
    check <- match.G2A(A = A[3:1, 3:1], G = dummymat, ord = TRUE)
  )
  check <- match.G2A(A = A[1:2, 1:2], G = dummymat, clean = TRUE, ord = FALSE, mism = TRUE,
                   RMdiff = FALSE)
  expect_equal(check$mismG, c("3", "4"))
})
test_that("traps work", {
  expect_error(
    check <- match.G2A(A = A, G = dummymat, clean = FALSE, ord = TRUE, mism = TRUE,
                       RMdiff = TRUE)
  )
  expect_error(
    check <- match.G2A(A = as.data.frame(A), G = dummymat)
  )
  expect_error(
    check <- match.G2A(A = A, G = as.data.frame(dummymat))
  )
  Awr <- A
  colnames(Awr) <- NULL
  expect_error(
  check <- match.G2A(A = Awr, G = dummymat)
  )
  Awr <- A
  rownames(Awr) <- NULL
  expect_error(
  check <- match.G2A(A = Awr, G = dummymat)
  )
  dummymatwr <- dummymat
  colnames(dummymatwr) <- NULL
  expect_error(
  check <- match.G2A(A = A, G = dummymatwr)
  )
  dummymatwr <- dummymat
  rownames(dummymatwr) <- NULL
  expect_error(
  check <- match.G2A(A = A, G = dummymatwr)
  )
  dummymatwr <- dummymat
  rownames(dummymatwr) <- 4:1
  expect_error(
  check <- match.G2A(A = A, G = dummymatwr)
  )
  Awr <- A
  rownames(Awr) <- 4:1
  expect_error(
  check <- match.G2A(A = Awr, G = dummymat)
  )
})
