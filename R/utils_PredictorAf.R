# this are subfunctions needed for teh computation. Some of them are taken from the 
# GitHub Repo of the iml package.

# used to compute the progress when AdaCB infill criterion is used
getProgressAdaCB = function(res.mbo, iter) {
  if (res.mbo$control$infill.crit$id == "adacb") {
    
    opdf = as.data.frame(res.mbo$opt.path)
    lambda.start = res.mbo$control$infill.crit$params$cb.lambda.start
    lambda.end = res.mbo$control$infill.crit$params$cb.lambda.end
    lambda = opdf[which(opdf$dob == iter), "lambda"]
    progress = (lambda - lambda.start) / (lambda.end - lambda.start)
    
  } else {
    progress = NULL
  }
}


# In order for the PredictorAf object to be created we need some internal function of the iml 
# package. All credits of such functions are reserved to iml authors. 
# Dowloaded from the iml Repo on Github.
# iml: 
#' @importFrom data.table data.table
Data <- R6::R6Class("Data",
                    public = list(
                      X = NULL,
                      y = NULL,
                      y.names = NULL,
                      feature.types = NULL,
                      feature.names = NULL,
                      n.features = NULL,
                      n.rows = NULL,
                      prob = NULL,
                      # Removes additional columns, stops when some are missing
                      match_cols = function(newdata) {
                        colnames_new <- colnames(newdata)
                        missing_columns <- setdiff(self$feature.names, colnames_new)
                        if (length(missing_columns) > 0) {
                          stop(sprintf("Missing columns: %s", paste(missing_columns, collapse = ", ")))
                        }
                        additional_columns <- setdiff(colnames_new, self$feature.names)
                        if (length(additional_columns) > 0) {
                          warning(sprintf("Dropping additional columns: %s", paste(additional_columns, collapse = ", ")))
                        }
                        newdata[self$feature.names]
                      },
                      sample = function(n = 100, replace = TRUE, prob = NULL, get.y = FALSE) {
                        if (is.null(prob) & !is.null(self$prob)) {
                          prob <- self$prob
                        }
                        indices <- sample.int(self$n.rows,
                                              size = n,
                                              replace = replace, prob = prob
                        )
                        if (get.y) {
                          cbind(self$X[indices, ], self$y[indices, ])
                        } else {
                          self$X[indices, ]
                        }
                      },
                      get.x = function(...) {
                        self$X
                      },
                      get.xy = function(...) {
                        cbind(self$X, self$y)
                      },
                      print = function() {
                        cat("Sampling from data.frame with", nrow(self$X), "rows and", ncol(self$X), "columns.")
                      },
                      initialize = function(X, y = NULL, prob = NULL) {
                        
                        assertDataFrame(X, all.missing = FALSE)
                        assertNamed(X)
                        if (length(y) == 1 & is.character(y)) {
                          assert_true(y %in% names(X))
                          self$y <- X[, y, drop = FALSE]
                          self$y.names <- y
                        } else if (inherits(y, "data.frame")) {
                          assertDataFrame(y, all.missing = FALSE, null.ok = TRUE, nrows = nrow(X))
                          self$y <- y
                          self$y.names <- colnames(self$y)
                          if (length(intersect(colnames(self$y), colnames(X))) != 0) {
                            stop("colnames of y and X have to be different.")
                          }
                        } else if (is.vector(y) | is.factor(y)) {
                          assert_vector(y, any.missing = FALSE, null.ok = TRUE, len = nrow(X))
                          self$y <- data.frame(.y = y)
                          self$y.names <- colnames(self$y)
                        }
                        self$X <- data.table::data.table(X[, setdiff(colnames(X), self$y.names), drop = FALSE])
                        if (ncol(self$X) == 1) stop("Only 1 feature was provided. The iml package is only useful and works for multiple features.")
                        self$prob <- prob
                        self$feature.types <- get.feature.type(unlist(lapply(self$X, class)))
                        self$feature.names <- colnames(self$X)
                        self$n.features <- ncol(self$X)
                        self$n.rows <- nrow(self$X)
                        names(self$feature.types) <- self$feature.names
                      }
                    )
)

# iml, used for Data
get.feature.type <- function(feature.class) {
  checkmate::assertCharacter(feature.class)
  
  feature.types <- c(
    "numeric" = "numerical",
    "integer" = "numerical",
    "character" = "categorical",
    "factor" = "categorical",
    "ordered" = "categorical"
  )
  
  stopifnot(all(feature.class %in% names(feature.types)))
  feature.types[feature.class]
}

# iml, used for Predictor and PredictorAdf
checkPrediction <- function(prediction, data) {
  checkmate::assert_data_frame(data)
  checkmate::assert_data_frame(prediction,
                               nrows = nrow(data), any.missing = FALSE,
                               types = c("numeric", "integerish", "factor")
  )
}

# iml, used for Predictor and PredictorAf. We actually need only a part of it, since PredictorAf only works for WrappedModel class
inferTaskFromModel <- function(model) {
  UseMethod("inferTaskFromModel")
}
inferTaskFromModel.WrappedModel <- function(model) {
  if (!requireNamespace("mlr")) {
    stop("Please install the mlr package.")
  }
  if (inherits(model, "WrappedModel")) {
    tsk <- mlr::getTaskType(model)
  }
  if (tsk == "classif") {
    if (model$learner$predict.type != "prob") {
      warning("Output seems to be class instead of probabilities. 
               Automatically transformed to 0 and 1 probabilities.
               You might want to set predict.type = 'prob' for Learner!")
    }
    return("classification")
  } else if (tsk == "regr") {
    return("regression")
  } else {
    stop(sprintf("mlr task type <%s> not supported", tsk))
  }
}

#iml, used for Predictor and PredictorAf
inferTaskFromPrediction <- function(prediction) {
  assert_true(any(class(prediction) %in%
                    c("integer", "numeric", "data.frame", "matrix", "factor", "character")))
  if (inherits(prediction, c("data.frame", "matrix")) && dim(prediction)[2] > 1) {
    "classification"
  } else if (inherits(prediction, c("factor", "character"))) {
    "classification"
  } else {
    "regression"
  }
}
