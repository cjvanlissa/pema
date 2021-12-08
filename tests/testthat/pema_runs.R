library(metafor)
df <- dat.baskerville2012
df$vi <- df$se^2
df[["se"]] <- NULL
df$yi <- df$smd
df[["smd"]] <- NULL
mods <- c("year", "score", "fumonths", "retention", "outcomes", "duration", "meetings", "hours")
df <- df[c("yi", "vi", mods)]
df <- missRanger::missRanger(df)

res <- tryCatch({
  pema::brma(yi ~. , df)
  }, error = function(e){NULL})
test_that("pema runs", {
  expect_true(!is.null(res))
})
