#' Percentage of inconsistent examples
#'
#' Compute the percentage of inconsistent examples, i.e., a set of examples
#' with the same attribute values and different class labels.
#'
#' @param data A data frame with discrete variables.
#' @param class.idx An optional numerical index for the position of the class
#' variable in the data frame. It is assumed the last variable in the data.
#'
#' @return A list with the percentage of inconsistent examples (pct) and a
#' nested list with the indexes of such examples (idx).
#' @export
#'
#' @examples
#' model <- NCAR(data = qb, noise = 0.1)
#' df <- generateDF(model$parameters, 100)
#' df <- df |> dplyr::mutate(Error = NULL,
#'   Class = ObservedClass, ObservedClass = NULL)
#' result <- inconsistency(df)
#' result$pct
#' result$idx
inconsistency <- function(data, class.idx = ncol(data)){

  if(all(c("Class", "ObservedClass", "Error") %in% names(data))){
    warning("The 'Error' and 'Class' variables are avoided to calculate the percentage of inconsistent examples.")
    data[,"Error"] <- NULL
    data[,"Class"] <- NULL
    class.idx = ncol(data)
  }
  if(!is.data.frame(data)){
    stop("Data must be a data frame.")
  }
  if(nrow(data) != nrow(Filter(is.factor, data))){
    stop("The dataset contains columns that are not factors.")
  }

  names(data)[class.idx] <- "Class"

  attribs <- names(data)
  attribs <- attribs[1:(length(attribs)-1)]

  row.id <- NULL
  set <- data |> dplyr::mutate(row.id = seq(1, nrow(data))) |>
    dplyr::group_by_at(attribs) |>
    dplyr::summarise(idx = list(row.id), .groups = "keep")

  inconsistency.idx = c()

  for (i in set$idx){
    if (length(unique(data[i,]$Class)) > 1){
      inconsistency.idx <- c(inconsistency.idx, i)
    }
  }

  list(pct = length(inconsistency.idx)/nrow(data), idx = list(inconsistency.idx))
}

#' Volume of the overlapping region
#'
#' Compute the ratio of overlapped examples in a discrete data frame.
#'
#' @param data A data frame with discrete variables.
#' @param class.idx An optional numerical index for the position of the class
#' variable in the data frame. It is assumed the last variable in the data.
#'
#' @return A float number representing the volume of the overlapping region.
#' @export
#'
#' @examples
#' data(qb)
#' F2 <- discF2(qb)
#' F2
discF2 <- function(data, class.idx = ncol(data)){

  if(!is.data.frame(data)){
    stop("Data must be a data frame.")
  }
  if(nrow(data) != nrow(Filter(is.factor, data))){
    stop("The dataset contains columns that are not factors.")
  }

  names(data)[class.idx] <- "Class"
  attributes <- colnames(data)[ !colnames(data) == 'Class']

  r = 1
  for(attribute in attributes){
    overlap = 0
    range = unique(data[[attribute]])
    for (value in range){
      cn <- length(unique(data[data[[attribute]]==value,"Class"]))
      if(cn > 1){
        overlap = overlap + 1
      }
    }
    r <- r*(overlap/length(range))
  }

  r
}

#' Collective Feature Efficiency
#'
#' Compute the ratio of examples with overlapping in a nominal data frame.
#'
#' @param data A data frame with discrete variables.
#' @param class.idx An optional numerical index for the position of the class
#' variable in the data frame. It is assumed the last variable in the data.
#'
#' @return A list with the value of F4 (value), a nested list of the best
#' attributes to separate examples (attribs), a nested list of indexes of
#' examples with overlapping (overlapping) and a nested list of indexes of
#' examples with no overlapping (nonOverlapping).
#' @export
#'
#' @examples
#' data(qb)
#' F4 <- discF4(qb[1:15,])
#' F4$value
#' F4$attribs
#' F4$overlapping
#' F4$nonOverlapping
discF4 <- function(data, class.idx = ncol(data)){

  if(!is.data.frame(data)){
    stop("Data must be a data frame.")
  }
  if(nrow(data) != nrow(Filter(is.factor, data))){
    stop("The dataset contains columns that are not factors.")
  }

  names(data)[class.idx] <- "Class"

  rows <- nrow(data)
  variables <- length(data)
  rownames(data) <- seq(rows)
  bestAttribs <- c()
  noOverlap <- c()

  while (variables > 1){

    attribs <- colnames(data)[ !colnames(data) == 'Class']

    bestSet <- c()
    bestAttrib <- ""

    for(attrib in attribs){
      range = unique(data[[attrib]])
      nonoverlapping <- c()

      for (value in range){
        split <- data[data[[attrib]]==value,]
        cn <- length(unique(split[,"Class"]))
        if(cn <= 1){
          nonoverlapping <- c(nonoverlapping, as.integer(rownames(split)))
        }
      }

      if(length(bestSet) < length(nonoverlapping)){
        bestSet <- nonoverlapping
        bestAttrib <- attrib
      }
    }

    if(bestAttrib != ""){
      bestAttribs <- c(bestAttribs, bestAttrib)
      noOverlap <- c(noOverlap, bestSet)
      data <- data[!(rownames(data) %in% bestSet), ]
      data[, bestAttrib] <- NULL
    }

    if(length(bestSet) == 0){
      break
    }

    variables <- length(data)
  }

  list(value = nrow(data)/rows, attribs = bestAttribs, overlapping = row.names(data), nonOverlapping = noOverlap)
}

#' Noise level
#'
#' Compute the noise level of data frame generated with the funtions of
#' noise models.
#'
#' @param data A data frame.
#'
#' @return A float number.
#' @export
#'
#' @examples
#' data(qb)
#' model <- NCAR(data = qb, noise = 0.4)
#' df <- generateDF(model$parameters, 1000)
#' result <- noiseLVL(df)
#' result
noiseLVL <- function(data){

  if(!is.data.frame(data)){
    stop("Data must be a data frame.")
  }
  if(!all(c("Error", "Class", "ObservedClass") %in% names(data))){
    stop("The dataset does not contain the columns to estimate the noise level.")
  }

  r <- data[data$Error == "True",] |> nrow()
  r <- (r * 100) / nrow(data)

  r
}
