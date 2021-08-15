source("R/utils_PredictorAf.R")
#' @title Predictor object that contains the Acquisition Function of a MBO process
#'
#' @description
#' `PredictorAf` has not been used in this project, but it is still provided as an additional resource
#' for users that want to assess the overall desirability of parameters with Acquisition functions other 
#' than the CB, for instance the EI. Tests for this method can be found at 
#' \url{https://github.com/croppi-f/iMBO/tree/master/tests/testthat}. `PredictorAf` is a extension / modification of iml package `Predictor`(inherits its behavior). It holds 
#' a `mlr` model (the surrogate model), the AF of a BO process and the candidate points to analyze.  
#' So far, PredictorAf works for the seven built in SingleObjective AF of mlrMBO.
#' PredictorAf is built on iml version ‘0.10.0’.Due to inherit, if latest versions of iml are released PredictorAf might generate errors.
#' 
#' @details
#' -`PredictorAf` has been created for regression tasks. For classification tasks it is not tested.
#' - `PredictorAf` class is created in a more general way, in order to be eventually compatible with other iml tools.
#' -Please use it only in combination with the `mlrMBO` package.
#' 
#' @name PredictorAf
#' @export
PredictorAf = R6::R6Class("PredictorAf",
  inherit = iml::Predictor,
  public = list(
    #' @description Create a PredictorAf object
    #' @param model\cr
    #' A machine learning model from `mlr` package, used as surrogate model in the 
    #' MBO process, and needed in order to measure the AF. 
    #' @param data [data.frame]\cr
    #' The data for which we want to measure the Af.
    #' @param y `character(1)` | [numeric] | [factor]\cr 
    #' The target vector or (preferably) the name of the target column (the id of the infill criterion) in the `data` argument.
    #' PredictorAf tries to infer the target automatically if not given. 
    #' @param batch.size `numeric(1)`\cr
    #' The maximum number of rows to be input the model for prediction at once.
    #' @param res.mbo [MBOSingleObjResult]\cr
    #' The result object of the MBO process.
    #' @param design [data.frame]\cr
    #' The design or training set used to fit the surrogate model. The argument must include also a column with the target value.
    #' Note, that this argument corresponds to the data argument of the Predictor Obejct in `iml`.
    #' @param iter `numeric(1)`\cr
    #' The iteration of interest within the MBO process
    initialize = function(
      # some arguments of the iml::Predictor 
      model = NULL, data = NULL, y = NULL, batch.size = 1000,
      # arguments needed for the Predictor of the Af
      res.mbo = NULL, design = NULL, iter = NULL
    ) {
      checkmate::assertClass(model, classes = "WrappedModel", null.ok = FALSE)
      checkmate::assertClass(res.mbo, classes = c("MBOSingleObjResult", "MBOResult"))
      opdf =  as.data.frame(res.mbo$opt.path)
      itmax = opdf[nrow(opdf), "dob"]
      checkmate::assertNumber(iter, lower = 1, upper = itmax)
      # extract useful informations
      par.set.mbo = res.mbo$opt.path$par.set
      pars.mbo = getParamIds(par.set.mbo, repeated = TRUE, with.nr = TRUE)
      y.name.mbo = res.mbo$control$y.name
      infill.mbo = res.mbo$control$infill.crit$id
      initial.design = res.mbo$final.opt.state$opt.problem$design
      acq.fun = res.mbo$control$infill.crit$fun
      
      # the design is used to compute the AF 
      if (is.null(design)) {
        stop("Please provide the design set")
      } else {
        # assertions works only for Single-Objective & Single-Point Proposals
        checkmate::assertDataFrame(design,
          any.missing = FALSE, 
          all.missing = FALSE,
          ncols = length(pars.mbo) + 1,
          nrows = nrow(initial.design) + iter - 1
        )
        checkmate::assertNames(names(design), identical.to = c(pars.mbo, y.name.mbo))
      }
      # y in this case is the id of the infill criterion and not the id of the target variable
      if (is.null(y)) {
        y = infill.mbo
        # y not always needed, so ignore when not in data (e.g. in the case of FEATURE EFFECT NOT NEEDED)
        if (is.character(y) && !(y %in% names(data))) {
          y = NULL
        }
      }
      # data used are the candidate points
      if (is.null(data)) {
        # notice that it extracts the data, but the value for the Af for those points is missing
        if (!is.null(res.mbo$seen.points)) data = res.mbo$seen.points[[iter]]
        else stop("Can't extract the candidate points. Provide them via data= or make sure to run the MBO with Savepts opt. methods")
      }
      checkmate::assertDataFrame(data, any.missing = FALSE, all.missing = FALSE)
      
      # initialize public fields
      self$data = Data$new(data, y = y)
      self$model = list(model)
      self$task = inferTaskFromModel(model)
      self$batch.size = batch.size
      self$res.mbo = res.mbo
      self$par.set = par.set.mbo
      self$control = res.mbo$control
      self$design = list(design)
      self$iter = iter
      self$progress = getProgressAdaCB(res.mbo, iter)
      # the prediction function is extracted from res.mbo, do not need an argument for that
      self$prediction.function = acq.fun
    },
    # taken form iml::Predictor and adjusted
    #' @description Predict the value of the AF for newdata, e.g. other candidate points.
    #' @template newdata
    predict = function(newdata) {
      checkmate::assert_data_frame(newdata)
      # Makes sure it's not a data.table
      newdata = as.data.frame(newdata)
      # make sure only features are used
      newdata = newdata[, intersect(
        self$data$feature.names,
        colnames(newdata)
      ), drop = FALSE]
      
      .prediction = self$prediction.function(
        points = newdata,
        models = self$model,
        control = self$control,
        par.set = self$par.set,
        designs = self$design,
        iter = self$iter,
        progress = self$progress,
        attributes = FALSE
      )
      .prediction = data.frame(.prediction)
      
      if (!private$predictionChecked) {
        checkPrediction(.prediction, newdata)
        private$predictionChecked = TRUE
      }
      # If S3 or function, we can only infer the task
      # once we see the predictions
      if (is.null(self$task)) {
        self$task = inferTaskFromPrediction(.prediction)
      }
      rownames(.prediction) = NULL
      data.frame(.prediction, check.names = FALSE)
    },
    #' @description Print the PredictorAf object.
    print = function() {
      cat(
        "Prediction task of surrogate model:", self$task,
        ", Infill criteria:", self$res.mbo$control$infill.crit$id, "\n"
      )
    },
    
    #' @field data [data.frame]\cr
    data = NULL,
    
    #' @field model (any)\cr
    model = NULL,
    
    #' @field batch.size `numeric(1)`\cr
    #' The number of rows to be input the model for prediction at once.
    batch.size = NULL,
   
    #' @field task `character(1)`\cr
    #'   The inferred prediction task: `"classification"` or `"regression"`. So far,
    #'   you will encounter only `"regression"`.
    task = NULL,
   
    # begin edited:
    #' @field res.mbo [MBOSingleObjResult]
    #' the results of the MBO process
    res.mbo = NULL,
   
    #' @field prediction.function
    #' The acquisition function extracted from the MBO process. So far, the seven
    #' built in AF for SinglObejctive are implemented.
    prediction.function = NULL,
   
    #' @field par.set [ParamSet]\cr
    #' The parameter set of the MBO process, describing different aspects of the
    #' objective function parameters
    par.set = NULL,
   
    #' @field control [MBOControl]\cr
    #' The control object of the MBO process
    control = NULL,
   
    #' @field design [data.frame]\cr
    #' The design set  used to compute the AF. 
    design = NULL,
   
    #' @field  iter `integer(1)`\cr
    #' The iteration of interest, in which we want to analyse the AF. Please make sure that the 
    #' corresponding surrogate model has been stored troughout the process (see store.model.at
    #' in makeMBOControl)
    iter = NULL,
   
    #' @field progress `numeric(1)`\cr
    #' A value between 0 and 1 indicating the progress of the optimization. Only needed in case of AdaCB,
    #' or other custom Adaptive Infill Criteria.
    progress = NULL
  ),
  private = list(
    predictionChecked = FALSE
  )
)
