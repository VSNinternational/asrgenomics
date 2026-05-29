datas <- pheno.apple |> head()
False <- "FALSE"
INDIV <- "INDIV"
INDIVIDUO <- "INDIVIDUO"
test_that("traps work", {
  expect_error(
    check.data_("datas", class_ = "matrix")
  )
  expect_error(
    check.data.mode_("datas", mode_ = "matrix")
  )
  expect_error(
    check.args_(
      data_ = datas, arg_ = INDIVIDUO, class_ = "integer")
  )
  expect_error(
    check.logical_("False")
  )
  expect_message(
    check.args_(
      data_ = datas, mandatory_ = F, arg_ = INDIV,
      class_ = "integer", class.action_ = "message", message_ = TRUE)
  )
})
