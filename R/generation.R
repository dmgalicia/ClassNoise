#' Generation of synthetic data
#'
#' Generate a data frame with an specific number of examples according the
#' distribution of a noise model.
#'
#' @param data A data frame with discrete variables.
#' @param size An integer representing the number of generated examples.
#'
#' @return A data frame.
#' @export
#'
#' @examples
#' data(qb)
#' model <- NCAR(data = qb, noise = 0.4S)
#' df <- generateDF(model$parameters, 10)
#' df
generateDF <- function(model, size){
  df <- bnlearn::rbn(model, size)
  df
}
