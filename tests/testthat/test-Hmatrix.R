dummymat <-
  matrix(c(1,0,0,
           0,1,.25,
           0,.25,1),
         byrow = TRUE, nrow = 3,
         dimnames = list(1:3, 1:3))
A <- diag(nrow = nrow(dummymat))
A[lower.tri(A)] <- A[upper.tri(A)] <- c(.5,.25,0)
rownames(A) <- colnames(A) <- rownames(dummymat)
Ginv <- G.inverse(G = dummymat)$Ginv
test_that("H matrix can be calculated", {
  H <- H.matrix(A = A, Ginv = Ginv, lambda = 0.90, sparseform = FALSE)
  expect_false(identical(H, dummymat))
})
