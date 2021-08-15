---
output:
  pdf_document: default
  html_document: default
  word_document: default
---
# ShapleyMBO
This repository contains all the relevant materials of my master's thesis project **Explaining Sequential Model Based Optimization** and is subdivided into following folders **R**, **dissertation** and **tests**.

1. **R**: contains all the functions used. The most important and "exported" functions are `ShapleyMBO`, `ShapleyMBO_mclapply`, `plotShapleyMBO` and `checkSampleSize`. This are ready to use for other users.

    - `ShapleyMBO`: does the computation and is the heart of the project. `ShapleyMBO_mclapply` is the parallel version thereof. This functions internally call other necessary **subfunctions** to compute the Shapley Values:
    - `ShapleyAf` and `PredictorAf` (see also `utils_PredictorAf.R`) are called to create the Shapley object when               `contribution = FALSE`.
    - `ShapleyMean` and `ShapleySe` are called separately when `contribution = TRUE`, i.e. when the LCB contributions are decomposed (topic of the master's thesis).
    - other utility functions are stored in `utils_ShapleyMBO.R`:
      - `getShapleyRes` is used to extract the results data frame from the Shapley object.
      - `mergeShapleyRes` is used to merge the results of Shapley Mean and Shapley Se and to compute than the actual and average cb.
      - `computePhiCb` is a subfunction of `mergeShapleyRes` and is used to compute the Shapley Value for the  the LCB.
    - `plotShapleyMBO`: displays the results of `ShapleyMBO`. This function internally calls other necessary **subfunctions**:
      - `plotShapleyCF` where CF stands for contribution FALSE and `plotShapleyCT` where CT stands for contribution TRUE. Each R script has other internal subfunctions, `<foo>.bar` and `<foo>.line`, used to display single iterations or so called desirbaility paths.
    - `checkSampleSize` is a method proposed to find a sufficiently high sample size argument for the SV computation.
    - other folders contain functions that were used for other specfic purposes:
      - **custom_function_hyp_ell**: `createMBOrunHypEll4d_noisy` is used to create the MBO runs specfic for the Hyper-Ellipsoid analysis and `plotShapleyHypEllmultiRun.bar.sd` is used to visualize ShapleyMBO results averaged over multiple BO runs using std.dev. as error bounds.
      - **Shapley_env**: contains similar `Shapley<foo>` functions, but the Shapley objects are stored as environments and not modified (standard `iml::Shapley()` object). This functions are only used to test `ShapleyMBO`.


2. **dissertation**: all data, code and figures used in the thesis.

    - **analysis**: has two main folders **hyper_ellipsoid** and **mlp** that contain code and data for the MBO runs and Shapley computation as well as the results of the analysis.
    - **figures**: figures used in the thesis.
    - **other_code**: some additional examples not part of the analysis. Each file is named with the respective section of the thesis.
    
3. **tests**: tests for some functions. Most important is script `test_ShapleyMBO.R`.
