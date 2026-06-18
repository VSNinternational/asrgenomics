dummymat <-
  matrix(c(1,.99,.99,
           .99,.99,.99,
           .99,.99,.99),
         byrow = TRUE, nrow = 3,
         dimnames = list(1:3, 1:3))
A <- diag(nrow = nrow(dummymat))
A[lower.tri(A)] <- A[upper.tri(A)] <- c(.5,.25,0)
rownames(A) <- colnames(A) <- rownames(dummymat)
suppressWarnings(Ginv <- G.inverse(G = dummymat, bend = T)$Ginv)
suppressWarnings(Ainv <- G.inverse(G = A)$Ginv)
test_that("H matrix can be calculated", {
  H <- H.inverse(A = A, Ginv = Ginv, lambda = 0.90, sparseform = TRUE)
  expect_equal(ncol(H), 3)
  H <- H.inverse(A = A, Ginv = Ginv, lambda = 0.90, inverse = F, sparseform = TRUE)
  expect_equal(ncol(H), 3)
})
test_that("traps work", {
  Awr <- A
  rownames(Awr) <- colnames(Awr) <- 3:5
  expect_error(
    H.inverse(A = Awr, Ginv = Ginv, sparseform = TRUE)
  )
  Awr <- A
  rownames(Awr) <- NULL
  expect_error(
    H.inverse(A = Awr, Ginv = Ginv, sparseform = TRUE)
  )
  Awr <- A
  colnames(Awr) <- NULL
  expect_error(
    H.inverse(A = Awr, Ginv = Ginv, sparseform = TRUE)
  )
  Awr <- A
  colnames(Awr)[1] <- 'nil'
  expect_error(
    H.inverse(A = Awr, Ginv = Ginv, sparseform = TRUE)
  )
  Ginvwr <- Ginv
  colnames(Ginvwr) <- NULL
  expect_error(
    H.inverse(A = A, Ginv = Ginvwr, sparseform = TRUE)
  )
  Ginvwr <- Ginv
  rownames(Ginvwr) <- NULL
  expect_error(
    H.inverse(A = A, Ginv = Ginvwr, sparseform = TRUE)
  )
  Ginvwr <- Ginv
  rownames(Ginvwr)[1] <- 'nil'
  expect_error(
    H.inverse(A = A, Ginv = Ginvwr, sparseform = TRUE)
  )
  Ainvwr <- Ainv
  rownames(Ainvwr) <- NULL
  expect_error(
    H.inverse(Ainv = Ainvwr, Ginv = Ginv, sparseform = TRUE)
  )
  Ainvwr <- Ainv
  colnames(Ainvwr) <- NULL
  expect_error(
    H.inverse(Ainv = Ainvwr, Ginv = Ginv, sparseform = TRUE)
  )
  Ainvwr <- Ainv
  colnames(Ainvwr)[1] <- 'nil'
  expect_error(
    H.inverse(Ainv = Ainvwr, Ginv = Ginv, sparseform = TRUE)
  )
  expect_warning(
    H.inverse(Ainv = Ainv, Ginv = Ginv, sparseform = TRUE)
  )
  Ainvwr <- Ainv
  diag(Ainvwr) <- 4
  expect_warning(
    H.inverse(A = Ainvwr, Ginv = Ginv, sparseform = TRUE)
  )
  expect_error(
    H.inverse(A = as.data.frame(A), Ginv = Ginv, sparseform = TRUE)
  )
  expect_error(
    H.inverse(A = A, Ginv = as.data.frame(Ginv), sparseform = TRUE)
  )
  expect_error(
    H.inverse(Ainv = as.data.frame(Ainv), Ginv = Ginv[3:1, 3:1])
  )
  expect_error(
    H.inverse(Ginv = Ginv)
  )
  expect_error(
    H.inverse(A = A, Ginv = Ginv, tau = 0, omega = 0)
  )
  expect_error(
    H.inverse(A = A, Ginv = Ginv, omega = -2)
  )
  expect_error(
    H.inverse(A = A, Ginv = Ginv, tau = -2)
  )
  expect_error(
    H.inverse(A = A, Ginv = Ginv, lambda = -2)
  )
})
