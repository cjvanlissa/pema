#' @title Compute model accuracy
#' @description Compute model accurary
#' @param fit A fit model object.
#' @param newdata New predictor data to evaluate the model fit, default: NULL.
#' @param observed Numeric vector of the dependent variable, default: NULL.
#' @param ymean Numeric. Mean to compare predicted values to, default: NULL.
#' @param olddata Original predictor data used to fit the model, default: NULL.
#' @param ... Additional arguments passed to 'predict' functions.
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ModelAccuracy
#' @keywords internal
#' @export
#' @importFrom methods hasArg
model_accuracy <-
  function(fit,
           newdata = NULL,
           observed = NULL,
           ymean = NULL,
           olddata = NULL,
           ...) {
    switch(
      class(fit)[1],
      brma = {
        if (is.null(newdata)) {
          predicted	<- pred_brma(fit)
        } else {
          predicted <- pred_brma(fit, data = newdata)
        }
      },
      rma.uni = {
        if (is.null(newdata)) {
          predicted	  <- predict(fit)$pred
          if (length(predicted) == 1) {
            predicted <- rep(predicted, length(observed))
          }
        } else {
          colnames(newdata) <- NULL
          predicted <- predict(fit, newmods = as.matrix(newdata))$pred
        }
      },

      MetaForest = {
        if (is.null(newdata)) {
          predicted <- fit$predictions
        } else {
          #colnames(newdata) <- NULL
          predicted <- predict(fit, data = newdata)$predictions
        }
      },
      FEmrt = {
        if (is.null(newdata)) {
          predicted	  <- suppressWarnings(predict(fit, newdata = fit$data)$g)
        } else {
          predicted <- suppressWarnings(predict(fit,
                                                newdata = newdata)$g)
        }
      },

      lma_it = {
        if (is.null(newdata)) {
          if(!hasArg("s")){
            predicted	<- predict(fit, newx = as.matrix(olddata), s= fit[[fit$use_lambda]], ...)
          } else {
            predicted	<- predict(fit, newx = as.matrix(olddata), ...)
          }
        } else {
          if(!hasArg("s")){
            predicted	<- predict(fit, newx = as.matrix(newdata), s= fit[[fit$use_lambda]], ...)
          } else {
            predicted	<- predict(fit, newx = as.matrix(newdata), ...)
          }
        }
      }
      ,
      cv.glmnet = {
        if (is.null(newdata)) {
          if(!hasArg("s")){
            predicted	<- predict(fit, newx = as.matrix(olddata), s= fit[["lambda.1se"]], ...)
          } else {
            predicted	<- predict(fit, newx = as.matrix(olddata), ...)
          }
        } else {
          if(!hasArg("s")){
            predicted	<- predict(fit, newx = as.matrix(newdata), s= fit[["lambda.1se"]], ...)
          } else {
            predicted	<- predict(fit, newx = as.matrix(newdata), ...)
          }
        }
      },

      {
        if (is.null(newdata)) {
          predicted	<- predict(fit, olddata, ...)
        } else {
          predicted	<- predict(fit, newdata, ...)
        }
      }
    )

    if (is.null(ymean)) ymean <- mean(observed)

    if (anyNA(predicted))
      message("Predictions for model of class ", class(fit)[1], " contained NAs. These were replaced with the value of 'ymean'.")
    predicted[is.na(predicted)] <- 0

    SS.total    <- sum((observed - ymean) ^ 2)
    SS.residual <- sum((observed - predicted) ^ 2)
    SS.regression <- sum((predicted - ymean) ^ 2)

    r.2	<- 1 - SS.residual / SS.total
    mse	<- SS.residual / length(observed)
    sd(predicted)
    if (sd(predicted) == 0) {
      r.actual.pred <- 0
    } else{
      r.actual.pred <- cor(observed, predicted)
    }

    c(r2 = r.2, mse = mse, r_actual_pred = r.actual.pred)
  }


pred_brma <- function(x, data = NULL, ...){
# Prepare data ------------------------------------------------------------
  if(is.null(data)){
    X <- fit$X
  } else {
    df <- data.frame(yi = 1, data) #I changed df <- data to df <- data.frame(yi=1, data)
    if(!is.null(x[["vi_column"]])){
      if(x[["vi_column"]] %in% names(df)) df[[x$vi_column]] <- NULL
    }
    if(!is.null(x[["study_column"]])){
      if(x[["study_column"]] %in% names(df)) df[[x[["study_column"]]]] <- NULL
    }
    mf <- call("model.frame",
               formula = x$formula,
               data = df,
               drop.unused.levels = TRUE)
    mf <- eval(mf, parent.frame())
    X <- mf[, -1, drop = FALSE]
  }

# Prepare coefs -----------------------------------------------------------

  coefs <- rstan::summary(x$fit)$summary[, "mean"]
  int <- coefs[names(coefs) == "Intercept"]
  coefs <- coefs[startsWith(names(coefs), "betas[")]

# Produce prediction ------------------------------------------------------
  int + rowSums(X * outer(rep.int(1L, nrow(X)), coefs))
}

