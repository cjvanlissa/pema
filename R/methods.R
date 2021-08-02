#' @method print brma
#' @export
print.brma <- function(x, ...){
  cl <- match.call()
  cl[["x"]] <- x[["model"]]
  eval.parent(cl)
}
