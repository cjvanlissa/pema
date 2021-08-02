#' @title rma for simulations
#' @description Tries rma() with several settings until it converges, or throws
#' an error if none of the settings result in a converging model. Used for
#' simulations only.
#' @param ... Arguments to be passed to rma()
#' @rdname rma_sim
#' @export
rma_for_sim <- function(...){
  #browser()
  args <- as.list(match.call()[-1])
  res <- try(do.call(rma, args), silent = TRUE)

  if(!is.atomic(res)) return(res)

  args$control <- list(stepadj=0.01, maxiter=1000)
  res <- try(do.call(rma, args), silent = TRUE)

  if(!is.atomic(res)) return(res)

  args$method <- "EB"
  res <- try(do.call(rma, args), silent = TRUE)

  if(!is.atomic(res)) return(res)

  args$method <- "ML"
  res <- try(do.call(rma, args), silent = TRUE)

  if(!is.atomic(res)) return(res)

  args$method <- "DL"
  res <- try(do.call(rma, args), silent = TRUE)

  if(!is.atomic(res)) return(res)

  data <- args$data
  largevalues <- which(abs(scale(data$yi)) > 3)
  data[largevalues, ] <- data[largevalues-1, ]
  args$data <- data
  args$control <- NULL
  res <- try(do.call(rma, args), silent = TRUE)

  if(is.atomic(res)){
    res <- NA
  }
  res
}
