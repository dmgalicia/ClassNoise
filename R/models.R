#' Noisy Completely At Random (NCAR) model
#'
#' This function generates and fits NCAR models with a specified noise level.
#' By default, it employs a hill-climbing algorithm to discover the Bayesian
#' Network structure. Alternatively, it can accept an \code{bn} object with
#' a pre-defined structure. The \code{data} argument can receive data frames
#' with numeric and categorical variables. The NCAR function also assumes the
#' last variable in the dataset serves as class and automatically renames it
#' as 'Class' during model creation.
#'
#' @param data A data frame containing the domain variables.
#' @param noise A number between 0 and 1 representing the noise level.
#' @param class.idx An optional numeric index for the position of the class
#' variable in the data frame. It is assumed the last variable in the data.
#' @param network A \code{bn} object representing the structure of the golden
#' Bayesian Network of the domain.
#'
#' @return A list with a \code{bn} object representing the structure and a
#' \code{bn.fit} object representing the parameterized model.
#' @export
#'
#' @examples
#' data(qb)
#'
#' # The following examples have P(Error = True) = 0.1
#' # Running the method without a predefined network
#' model <- NCAR(data = qb, noise = 0.1)
#' model$parameters$Error
#' # Running the method with a predefined network
#' network <- bnlearn::hc(qb)
#' fit <- bnlearn::bn.fit(network, qb)
#' model <- NCAR(data = qb, noise = 0.1, network = fit)
#' model$parameters$Error
NCAR <- function(data, noise, class.idx = ncol(data), network = NULL){

  if(!is.data.frame(data)){
    stop("Data must be a data frame.")
  }
  if(!is.factor(data[,class.idx])){
    stop("The class must be a factor.")
  }

  names(data)[class.idx] <- "Class"

  labels <- sort(unique(as.vector(data[, "Class"])))
  num.labels <- length(labels)
  attribs <- names(data)
  attribs <- attribs[attribs != "Class"]
  num.attribs <- length(attribs)
  struct <- NULL
  data[, "Error"] <- factor(rep(c("False", "True"), length.out = nrow(data)),levels = c("False", "True"))
  data[, "ObservedClass"] <- data[, "Class"]

  if(is.null(network)){

    bl <- rbind(data.frame(from = names(data), to = rep("ObservedClass", times = ncol(data))),
                data.frame(from = names(data), to = rep("Error", times = ncol(data))),
                data.frame(from = c("ObservedClass", "Error"), to = c("Error", "ObservedClass")),
                data.frame(from = rep("ObservedClass", times = ncol(data)), to = names(data)),
                data.frame(from = rep("Error", times = ncol(data)), to = names(data)))

    struct <- bnlearn::hc(data, restart = 100, perturb = 10, blacklist = bl)

    bn.arcs <- bnlearn::arcs(struct)
  }else{
    if(!all(bnlearn::nodes(network) %in% names(data))){
      stop("The network data frame does not contain the same variables as the provided data.")
    }

    struct <- bnlearn::empty.graph(nodes = c(bnlearn::nodes(network), "Error", "ObservedClass"))
    bn.arcs <- bnlearn::arcs(network)
  }

  noise.arcs <- data.frame(from = c("Class", "Error"), to = c("ObservedClass", "ObservedClass"))
  bnlearn::arcs(struct) <- rbind(bn.arcs, noise.arcs)

  param <- bnlearn::bn.fit(struct, data)

  param$Error = matrix(c(1-noise, noise), ncol = 2, dimnames = list(NULL, c("False", "True")))

  cpt.oc.F <- diag(num.labels)
  cpt.oc.T <- ((diag(num.labels) - 1) * (-1)) * (1 / (num.labels - 1))
  cpt.oc <- c(as.vector(cpt.oc.F), as.vector(cpt.oc.T))
  dim(cpt.oc) <- c(num.labels, num.labels, 2)
  oc.names <- list()
  oc.names[["ObservedClass"]] <- labels
  oc.names[["Class"]] <- labels
  oc.names[["Error"]] <- c("False","True")
  dimnames(cpt.oc) <- oc.names
  param$ObservedClass <- cpt.oc

  list(structure = struct, parameters = param)
}


#' Noise At Random (NAR) model
#'
#' This function generates and fits NAR models where each class label has a
#' specific noise level. Noise levels are introduced as a numeric vector
#' which follows the order of the class factor. By default, it employs a
#' hill-climbing algorithm to discover the Bayesian Network structure.
#' Alternatively, it can accept an \code{bn} object with a pre-defined
#' structure. The \code{data} argument can receive data frames with numeric
#' and categorical variables. The NAR function also assumes the last variable
#' in the dataset serves as class and automatically renames it as 'Class'
#' during model creation.
#'
#' @param data A data frame containing the domain variables.
#' @param noise A number or a numeric vector of \code{n} elements between 0
#' and 1 representing the percentage of noise affecting each label. If a
#' number is provided all clases have the same probability of error
#' @param class.idx An optional numeric index for the position of the class
#' variable in the data frame. It is assumed the last variable in the data.
#' @param network A \code{bn} object representing the structure of the golden
#' Bayesian Network of the domain.
#'
#' @return A list with a \code{bn} object representing the structure and a
#' \code{bn.fit} object representing the parameterized model.
#' @export
#'
#' @examples
#' data(qb)
#'
#' # The following examples have P(Error = True | Class = B) = 0.5
#' # and P(True | Error = Class = NB) = 0
#' # Running the method without a predefined network
#' model <- NAR(data = qb, noise = c(0.5, 0))
#' model$parameters$Error
#' # Running the method with a predefined network
#' network <- bnlearn::hc(qb)
#' fit <- bnlearn::bn.fit(network, qb)
#' model <- NAR(data = qb, noise = c(0.5, 0), network = fit)
#' model$parameters$Error
NAR <- function(data, noise, class.idx = ncol(data), network = NULL){

  if(!is.data.frame(data)){
    stop("Data must be a data frame.")
  }
  if(!is.factor(data[,class.idx])){
    stop("The class must be a factor.")
  }

  names(data)[class.idx] <- "Class"
  labels <- sort(unique(as.vector(data[, "Class"])))
  num.labels <- length(labels)
  attribs <- names(data)
  attribs <- attribs[attribs != "Class"]
  num.attribs <- length(attribs)
  struct <- NULL
  data[, "Error"] <- factor(rep(c("False", "True"), length.out = nrow(data)),levels = c("False", "True"))
  data[, "ObservedClass"] <- data[, "Class"]

  if(is.null(network)){

    bl <- rbind(data.frame(from = names(data), to = rep("ObservedClass", times = ncol(data))),
                data.frame(from = names(data), to = rep("Error", times = ncol(data))),
                data.frame(from = c("ObservedClass", "Error"), to = c("Error", "ObservedClass")),
                data.frame(from = rep("ObservedClass", times = ncol(data)), to = names(data)),
                data.frame(from = rep("Error", times = ncol(data)), to = names(data)))

    struct <- bnlearn::hc(data, restart = 100, perturb = 10, blacklist = bl)

    bn.arcs <- bnlearn::arcs(struct)
  }else{
    if(!all(bnlearn::nodes(network) %in% names(data))){
      stop("The network data frame does not contain the same variables as the provided data.")
    }

    struct <- bnlearn::empty.graph(nodes = c(bnlearn::nodes(network), "Error", "ObservedClass"))
    bn.arcs <-  bnlearn::arcs(network)
  }

  noise.arcs <- data.frame(from = c("Class", "Error", "Class"),
                           to = c("ObservedClass", "ObservedClass", "Error"))
  bnlearn::arcs(struct) <- rbind(bn.arcs, noise.arcs)

  param <- bnlearn::bn.fit(struct, data)

  error.names <- list()
  error.names[["Error"]] <- c("False","True")
  error.names[["Class"]] <- labels

  cpt.error <- c()
  if ( length(noise) == 1 ) {
    for (i in 1:num.labels) {
      cpt.error <- append(cpt.error, c(1-noise, noise))
    }
  } else if ( length(noise) == num.labels){
    for (i in 1:num.labels) {
      cpt.error <- append(cpt.error, c(1-noise[i], noise[i]))
    }
  } else {
    stop("Noise parameter must be a number or a numeric vector of n elements, where n is the number of class labels.")
  }
  cpt.error <- matrix(cpt.error, ncol = num.labels, dimnames = error.names)
  param$Error <- cpt.error

  cpt.oc.F <- diag(num.labels)
  cpt.oc.T <- ((diag(num.labels) - 1) * (-1)) * (1 / (num.labels - 1))
  cpt.oc <- c(as.vector(cpt.oc.F), as.vector(cpt.oc.T))
  dim(cpt.oc) <- c(num.labels, num.labels, 2)
  oc.names <- list()
  oc.names[["ObservedClass"]] <- labels
  oc.names[["Class"]] <- labels
  oc.names[["Error"]] <- c("False","True")
  dimnames(cpt.oc) <- oc.names
  param$ObservedClass <- cpt.oc

  list(structure = struct, parameters = param)
}


#' Noisy Not At Radom (NNAR) model
#'
#' This function generates and fits NNAR models where each instance of an
#' attribute set and the class has a specific noise level. Noise levels are
#' introduced as a numeric vector which follows the order of the Cartesian
#' product of their factors. The number of parameters increases with the arity
#' of the attribute set and the possible values of each attribute. By default,
#' it employs a hill-climbing algorithm to discover the Bayesian Network
#' structure within the domain. Alternatively, it can accept an \code{bn}
#' object with an pre-defined structure. The \code{data} argument just receives
#' categorical data frames. The NNAR function also assumes the last variable
#' in the dataset serves as class and automatically renames it as 'Class'
#' during model creation.
#'
#' @param data A data frame containing categorical variables.
#' @param attrib.set A vector of strings representing the attributes related
#' with Error.
#' @param noise A number or a numeric vector of \code{n} elements between 0
#' and 1 representing the percentage of noise affecting each label. If a number
#' is provided all possible scenarios have the same probability of error
#' @param class.idx An optional numeric index for the position of the class
#' variable in the data frame. It is assumed the last variable in the data.
#' @param network A \code{bn} object representing the structure of the golden
#' Bayesian Network of the domain.
#'
#' @return A list with a \code{bn} object representing the structure and a
#' \code{bn.fit} object representing the parameterized model.
#' @export
#'
#' @examples
#' data(qb)
#'
#' # This examples have P(Error = True | IR = P, Class = B) = 0.4,
#' # P(Error = True | IR = A, Class = B) = 0),
#' # P(Error = True | IR = N, Class = B) = 0.2,
#' # P(Error = True | IR = A, Class = NB) = 0),
#' # P(Error = True | IR = N, Class = NB) = 0.1,
#' # and P(Error = True | IR = A, Class = NB) = 0.3)
#' # Running the method without a predefined network
#' model <- NNAR(data = qb, attrib.set = c("IR"),
#'  noise = c(0.4, 0, 0.2, 0, 0.1, 0.3))
#' model$parameters$Error
#' # Running the method with a predefined network
#' network <- bnlearn::hc(qb)
#' fit <- bnlearn::bn.fit(network, qb)
#' model <- NNAR(data = qb, attrib.set = c("IR"),
#'   noise = c(0.4, 0, 0.2, 0, 0.1, 0.3), network = fit)
#' model$parameters$Error
NNAR <- function(data, attrib.set, noise, class.idx = ncol(data), network = NULL){

  if(!is.data.frame(data)){
    stop("Data must be a data frame.")
  }
  flag <- F
  for (i in attrib.set){
    if (is.numeric(data[,i])){flag <- T}
  }
  if(flag){
    stop("The attribute set cannot contain continuous attributes.")
  }

  names(data)[class.idx] <- "Class"

  attrib.set <- rev(attrib.set)
  labels <- sort(unique(as.vector(data[, "Class"])))
  num.labels <- length(labels)
  attribs <- names(data)
  attribs <- attribs[attribs != "Class"]
  num.attribs <- length(attribs)
  struct <- NULL
  data[, "Error"] <- factor(rep(c("False", "True"), length.out = nrow(data)),levels = c("False", "True"))
  data[, "ObservedClass"] <- data[, "Class"]

  if(is.null(network)){

    bl <- rbind(data.frame(from = names(data), to = rep("ObservedClass", times = ncol(data))),
                data.frame(from = names(data), to = rep("Error", times = ncol(data))),
                data.frame(from = c("ObservedClass", "Error"), to = c("Error", "ObservedClass")),
                data.frame(from = rep("ObservedClass", times = ncol(data)), to = names(data)),
                data.frame(from = rep("Error", times = ncol(data)), to = names(data)))

    struct <- bnlearn::hc(data, restart = 100, perturb = 10, blacklist = bl)

    bn.arcs <- bnlearn::arcs(struct)
  }else{
    if(!all(bnlearn::nodes(network) %in% names(data))){
      stop("The network data frame does not contain the same variables as the provided data.")
    }

    struct <- bnlearn::empty.graph(nodes = c(bnlearn::nodes(network), "Error", "ObservedClass"))
    bn.arcs <- bnlearn::arcs(network)
  }

  noise.arcs <- data.frame(from = c("Class", "Error", "Class"), to = c("ObservedClass", "ObservedClass", "Error"))
  for (i in attrib.set) { noise.arcs[nrow(noise.arcs) + 1,] <- c(i, "Error") }
  bnlearn::arcs(struct) <- rbind(bn.arcs, noise.arcs)

  param <- bnlearn::bn.fit(struct, data)

  error.names <- list()
  error.names[["Error"]] <- c("False","True")
  dimensions <- c(num.labels)
  for (i in attrib.set){
    error.names[[i]] <- levels(data[, i])
    dimensions <- c(dimensions, length(levels(data[, i])))
  }
  error.names[["Class"]] <- levels(data[, class.idx])
  dimensions <- c(dimensions, length(levels(data[, class.idx])))

  cpt.error <- c()
  if ( length(noise) == 1 ) {
    for (i in seq(prod(dimensions)/2)) {
      cpt.error <- append(cpt.error, c(1-noise, noise))
    }
  } else if ( length(noise) == prod(dimensions)/2){
    for (i in seq(prod(dimensions)/2)) {
      cpt.error <- append(cpt.error, c(1-noise[i], noise[i]))
    }
  } else {
    stop("Noise parameter must be a number or a numeric vector of n elements, where n is arity of the cartesian product of the values of attibutes related with Error and the class labels.")
  }
  cpt.error <- array(cpt.error, dim = dimensions, dimnames = error.names)
  param$Error <- cpt.error

  cpt.oc.F <- diag(num.labels)
  cpt.oc.T <- ((diag(num.labels) - 1) * (-1)) * (1 / (num.labels - 1))
  cpt.oc <- c(as.vector(cpt.oc.F), as.vector(cpt.oc.T))
  dim(cpt.oc)= c(num.labels, num.labels, 2)
  oc.names <- list()
  oc.names[["ObservedClass"]] <- labels
  oc.names[["Class"]] <- labels
  oc.names[["Error"]] <- c("False","True")
  dimnames(cpt.oc) <- oc.names
  param$ObservedClass <- cpt.oc

  list(structure = struct, parameters = param)
}
