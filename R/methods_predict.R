#' @method predict brma
#' @export
predict.brma <- function(object, newdata, type = c("mean", "median", "samples"), ...){
  if (missing(newdata) || is.null(newdata)) {
    X <- object$X
  } else {
    if(isTRUE(all(colnames(object$X) %in% colnames(newdata)))){
      X <- newdata[, colnames(object$X), drop = FALSE]
    } else {
      X <- try({
        form <- object$formula
        form[2] <- NULL
        if(!is.null(object[["vi_column"]])){
          if(object$vi_column %in% names(newdata)) newdata[[object$vi_column]] <- NULL
        }
        if(!is.null(object[["study_column"]])){
          if(object$study_column %in% names(newdata)) newdata[[object$study_column]] <- NULL
        }
        # Make model matrix
        mf <- call("model.frame")
        #mf <- mf[c(1L, match(c("formula", "subset", "na.action"), names(mf), 0L))]
        mf[["formula"]] <- form
        mf[["data"]] <- newdata
        mf$drop.unused.levels <- FALSE
        mf <- eval(mf, parent.frame())
        #X <- mf[, -1, drop = FALSE]
        #mt <- attr(mf, "terms")
        mt <- object$terms
        model.matrix(mt, mf)
      })
      if(inherits(X, "try-error")) stop("Could not construct a valid model.frame from argument 'newdata'.")
    }
  }
  switch(type[1],
         "mean" = .pred_brma_summary(object, X, parcol = "mean"),
         "median" = .pred_brma_summary(object, X, parcol = "50%"),
         "samples" = .pred_brma_samples(object, X))
}



.pred_brma_samples <- function(object, X, ...){
  sums <- summary(object$fit)$summary
  sim <- object$fit@sim
  keepthese <- c(which(sim$fnames_oi == "sd_1[1]"),
                 which(sim$fnames_oi == "Intercept"),
                 grep("^b\\[\\d+\\]$", sim$fnames_oi))

  sums <- sums[keepthese, , drop = FALSE]

  samps <- sapply(keepthese, .extract_samples, sim = sim)
  preds <- apply(samps[, -1], 1, function(thisrow){ rowSums(X %*% diag(thisrow)) })
  attr(preds, "tau") <- samps[,1]
  class(preds) <- c("brma_preds", class(preds))
  return(preds)
}

.pred_brma_summary <- function(object, X, parcol, ...){
  # Prepare coefs -----------------------------------------------------------
  coefs <- object$coefficients
  numpars <- sum(grepl("^b\\[\\d+\\]", object$fit@sim$fnames_oi))
  coefs <- coefs[1:(numpars+1), parcol]
  # Produce prediction ------------------------------------------------------
  unname(rowSums(X %*% diag(coefs)))
}
