#' Convert an object to stanfit
#'
#' Create a `stanfit` object from an object for which a method exists,
#' so that all methods for `stanfit` objects can be used.
#' @param x An object for which a method exists.
#' @param ... Arguments passed to or from other methods.
#' @return An object of class `stanfit`, as documented in [rstan::stan].
#' @export
#' @examples
#' stanfit <- "a"
#' class(stanfit) <- "stanfit"
#' brmaobject <- list(fit = stanfit)
#' class(brmaobject) <- "brma"
#' converted <- as.stan(brmaobject)
#' @importFrom stats rbinom rnorm rt
#' @importFrom sn rsn
as.stan <- function(x, ...){
  UseMethod("as.stan", x)
}

#' @method as.stan brma
#' @export
as.stan.brma <- function(x, ...){
  if(is.null(x[["fit"]])) stop("Could not coerce object to class 'stanfit'.")
  if(!inherits(x[["fit"]], "stanfit")) stop("Could not coerce object to class 'stanfit'.")
  out <- x[["fit"]]
  as_atts <- names(x)[!names(x) %in% c("fit")]
  for(thisatt in as_atts){
    attr(out, which = thisatt) <- x[[thisatt]]
  }
  attr(out, which = "type") <- "brma"
  return(out)
}
