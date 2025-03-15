
#' @title Test formulas for modeling tile accessibility (from MOCHA)
#'
#' @description \code{pilot_scATAC} Functions for testing out zero-inflated GLMM formulas on the data from MOCHA. Runs a given formula on a subset of the data, and returns the model results. This is meant to help during the model selection process. \code{\link[glmmTMB]{glmmTMB}}. 
#'
#' @param atacSE A MOCHA Tile-by-Sample Object (SummarizedExperiment)generated from getSampleTileMatrix within \code{\link[MOCHA]{MOCHA}}. 
#' @param cellPopulation Name of a cell type within the atacSE
#' @param modelFormula The formula for the continuous data that should be used within glmmTMB. It should be in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the atacSE metadata, except for CellType, FragNumber and CellCount, which will be extracted from the atacSE.
#'   modelFormula must start with 'exp' as the response.
#'   See \link[glmmTMB]{glmmTMB}.
#' @param ziFormula The formula for the zero-inflated data that should be used within glmmTMB. It should be in the
#'   format ( ~ factors). All factors must be found in column names
#'   of the atacSE colData metadata, except for CellType, FragNumber and CellCount, which will be extracted from the atacSE.
#'   FragNumber and CellCounts will be log10 normalized within the function. 
#' @param zi_threshold Zero-inflated threshold ( range = 0-1), representing the fraction of samples with zeros. At or above this threshold, the zero-inflated modeling kicks in.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores Number of cores to provide to glmmTMB for parallelization. 
#' @param featureSample integer or vector of characters describing the features to test. If integer, stratified sampling will occur to identify features to test, with a minimum of 10. If a character vector, then all strings must align to features you want to test. 
#'
#' @return results model results
#'
#'
#'
#' @examples
#' \dontrun{
#'   modelList <- pilot_scATAC(STM, 
#'                            'CD16 Mono',
#'                            exp~ Age + Sex + days_since_symptoms + (1|PTID),
#'                             ~ FragNumber, verbose = TRUE )
#' }
#'
#' @export
#' @keywords pilot_modeling

pilot_scATAC <- function(atacSE,
                        cellPopulation = NULL,
                        modelFormula = NULL,
                        ziFormula = ~0,
                        zi_threshold = 0,
                        verbose = FALSE,
                        numCores = 1,
                        featureSample = 10) {

  if(isChAIObject(atacSE, type = 'data', returnType = TRUE) != 'scATAC'){

    stop('atacSE is not a MOCHA SampleTileObject.')

  }


  if (!requireNamespace("MOCHA", quietly = TRUE)) {
      stop(
      "Package 'MOCHA' is required for pilot_scATAC. ",
      "Please install 'MOCHA' to proceed."
      )
  }
    
  if (!methods::is(modelFormula,'formula')) {
    stop("modelFormula was not provided as a formula.")
  }
    
    if (!methods::is(ziFormula,'formula')) {
    stop("ziFormula was not provided as a formula.")
  }
  if (zi_threshold < 0 | zi_threshold > 1 | ! is.numeric(zi_threshold)) {
    stop("zi_threshold must be between 0 and 1.")
  }

 if (length(cellPopulation) > 1) {
    stop(
      "More than one cell population was provided. ",
      "cellPopulation must be length 1. To run over multiple cell types, ",
      "user MOCHA to run combineSampleTileMatrix() to produce a new combined atacSE and set ",
      "cellPopulation = 'counts'."
    )
  } else if (
    (!cellPopulation %in% names(SummarizedExperiment::assays(atacSE)))
  ) {
    stop("cellPopulation was not found within atacSE.")
  } else if(cellPopulation == 'counts'){
    newObj <- atacSE
  }else{
    newObj <- MOCHA::combineSampleTileMatrix(MOCHA::subsetMOCHAObject(atacSE, subsetBy = 'celltype', groupList = cellPopulation, subsetPeaks = TRUE))
  }

  modelList <- .pilotModels_generic(SE_Object = newObj,
                        modelFormula = modelFormula,
                        ziFormula = ziFormula,
                        zi_threshold = zi_threshold,
                        family = stats::gaussian(),
                        modality = 'scATAC_Model',
                        verbose = verbose,
                        numCores = numCores,
                        featureSample = featureSample)
  
  return(modelList)
}


#' @title Test formulas for modeling gene expression
#'
#' @description \code{pilot_scRNA} Function for testing out GLM formulas on the pseudobulked RNA via ChAI.
#'   Runs a given formula on a subset of the data, and returns the model results. This is meant to help during the model selection process. \code{\link[glmmTMB]{glmmTMB}}. 
#'
#' @param rnaSE A SummarizedExperiment object generated from makePseudobulkRNA. 
#' @param cellPopulation Name of a cell type within the rnaSE
#' @param modelFormula The formula for the continuous data that should be used within glmmTMB. It should be in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the rnaSE metadata, except for CellType, FragNumber and CellCount, which will be extracted from the rnaSE.
#'   modelFormula must start with 'exp' as the response.
#'   See \link[glmmTMB]{glmmTMB}.
#' @param family String. Can be 'negativeBinomial1', 'negativeBinomial2', or 'poisson'. Default is "negativeBinomial2".
#' @param detectionThreshold A number between 0 and 1, representing the mean detection rate threshold for a given gene to be modeled. 
#'   This detection rate is calculated for each sample and cell type during makePseudobulkRNA, and represents the percentage of cells that have a transcript for a given gene. Over all samples, the average has to be above this to be modeled. Default is 0.01. 
#' @param expressionThreshold A number greater than zero, representing the expression threshold for modeling. A given gene, on average across all samples, expressed above this threshold. The default is 0. 
#' @param cellCountThreshold The minimum number of cells in a given pseudobulk for the pseudobulk to be included in analysis. If fewer than this number of cells are found, then the sample will be dicarded The number of cells within the pseudobulked scRNA. Default is 10 cells. 
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores Number of cores to provide to glmmTMB for parallelization. 
#' @param featureSample integer or vector of characters describing the features to test. If integer, stratified sampling will occur to identify features to test, with a minimum of 10. If a character vector, then all strings must align to features you want to test.  
#'
#' @return results model results
#'
#'
#'
#' @examples
#' \dontrun{
#'   modelList <- pilot_scRNA(rnaSE, 
#'                            'CD16 Mono',
#'                            exp~ Age + Sex + days_since_symptoms + (1|PTID),
#'                              verbose = TRUE )
#' }
#'
#' @export
#' @keywords pilot_modeling
#'

pilot_scRNA <- function(rnaSE,
                        cellPopulation = NULL,
                        modelFormula = NULL,
                        ziFormula = ~0,
                        zi_threshold = 0,
                        verbose = FALSE,
                        family = "negativeBinomial2",
                        detectionThreshold = 0.1,
                        expressionThreshold = 0,
                        cellCountThreshold = 0,
                        numCores = 1,
                        featureSample = 100) {

  if(isChAIObject(rnaSE, type = 'data', returnType = TRUE) != 'scRNA'){

    stop('rnaSE is not a ChAI scRNA Object (normalized, pseudobulked via ChAI.)')

  }

    
  if (!methods::is(modelFormula,'formula')) {
    stop("modelFormula was not provided as a formula.")
  }
    
    if (!methods::is(ziFormula,'formula')) {
    stop("ziFormula was not provided as a formula.")
  }

 if (length(cellPopulation) > 1) {
    stop(
      "More than one cell population was provided. ",
      "cellPopulation must be length 1. To run over multiple cell types, ",
      "run flattenChAI() to produce a new combined rnaSE and set ",
      "cellPopulation = 'counts'."
    )
  } else if (
    (!cellPopulation %in% names(SummarizedExperiment::assays(rnaSE)))
  ) {
    stop("cellPopulation was not found within rnaSE.")
  } else if(cellPopulation == 'counts'){
    newRNA <- rnaSE
  }else{
    #newObj <- flattenChAI(rnaSE)
    #Needs to be generated, analoguous to combineSampleTileMatrix

    newRNA = flattenChAI(rnaSE,  cellPopulations = cellPopulation, metadataT = TRUE)

  }

  if(tolower(family) == 'negativebinomial2'){ 
    family = glmmTMB::nbinom2()
  }else if(tolower(family)== 'negativebinomial1'){
    family = glmmTMB::nbinom1()
  }else if(tolower(family) == 'poisson'){
    family = stats::poisson()
  }else{
    stop('family not recognized.')
  }

  if(!is.numeric(detectionThreshold)){
    stop('detectionThreshold must be numeric.')
  }

  if(!is.numeric(expressionThreshold)){
    stop('expressionThreshold must be numeric.')
  }

  ### Filter by minimum detection rate to help speed up modeling. No point in modeling lowly expressed genes. 
  detectMean <- rowMeans(SummarizedExperiment::assays(rnaSE@metadata$detectionRate[rownames(newRNA),])[[cellPopulation]])
  expressMean <- rowMeans(SummarizedExperiment::assays(newRNA)[[1]])
  if(sum(detectMean > detectionThreshold & expressMean > expressionThreshold) == 0){
    stop('No genes pass the detectionThreshold and expressionThreshold. Please adjust thresholding.')
  }
    
  newRNA = newRNA[detectMean > detectionThreshold & expressMean > expressionThreshold,]

  if(sum(newRNA$CellCounts >=  cellCountThreshold) < 3){
  
      stop('Fewer than 3 samples passed the CellCounts threshold. Modeling not possible.')
      
  }else if(any(!newRNA$CellCounts >=  cellCountThreshold)){
  
      removeSamples = sum(newRNA$CellCounts <  cellCountThreshold)
      warning(stringr::str_interp('${removeSamples} Samples removed due to cell counts below cellCountThreshold'))
   
      newRNA <- newRNA[,newRNA$CellCounts >=  cellCountThreshold]
              
  }

  modelList <- .pilotModels_generic(SE_Object = newRNA,
                        modelFormula = modelFormula,
                        ziFormula = ziFormula,
                        family = family,
                        modality = 'scRNA_Model',
                        zi_threshold = zi_threshold,
                        verbose = verbose,
                        numCores = numCores,
                        featureSample = featureSample)
  
  return(modelList)
}




#' @title Test formulas for modeling ChromVAR motif Z-scores
#'
#' @description \code{pilot_chromVAR} Function for testing out GLM formulas on the ChromVAR Motif Z-scores generated via ChAI and MOCHA
#'   Runs a given formula on a subset of the data, and returns the model results. This is meant to help during the model selection process. \code{\link[glmmTMB]{glmmTMB}}. 
#'
#' @param MotifObj A SummarizedExperiment-type object generated from makeChromVAR.
#' @param cellPopulation Name of a cell type within the ChAI ChromVAR Object
#' @param modelFormula The formula to use with glmmTMB, in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the ExperimentObj metadata.
#' @param featureSample integer or vector of characters describing the features to test. 
#' If integer, stratified sampling will occur to identify features to test, with a minimum of 10. 
#' If a character vector, then all strings must align to features you want to test. 
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return modelList a list of outputs from glmmTMB
#'
#' @keywords pilot_modeling
#' @export
pilot_chromVAR <- function(MotifObj,
                      cellPopulation = NULL,
                      modelFormula,
                      numCores = 1,
                      featureSample = 100,
                      verbose = FALSE) {

  if(isChAIObject(chromSE, type = 'data', returnType = TRUE) != 'ChromVAR'){

    stop('chromSE is not a ChAI ChromVAR Object, created via makeChromVAR.')

  }
    
  if (length(cellPopulation) > 1) {
    stop(
      "More than one assay was provided. ",
      "assayName must be length 1. To run over multiple assays/cell types, ",
      "combine the matrices into one and wrap them in a SummarizedExperiment. Then set assayName to ",
      "the name of that matrix in the obejct."
    )
  } else if (
    !cellPopulation %in% names(SummarizedExperiment::assays(MotifObj))
  ) {
    stop("cellPopulation was not found within MotifObj provided.")
  }
    
  subObj <- flattenChAI(MotifObj, cellPopulations = cellPopulation) 

  modelList <- .pilotModels_generic(SE_Object = subObj,
                        modelFormula = modelFormula,
                        ziFormula = ~0,
                        zi_threshold = 0,
                        family = stats::gaussian(),
                        modality = 'ChromVAR_Model',
                        verbose = verbose,
                        numCores = numCores,
                        featureSample = featureSample)
  return(modelList)

}






#' @title Test formulas for modeling a general modality
#'
#' @description \code{pilot_General} Generalized Function for testing out zero-inflated and standard GLM formula on a general modality imported via \code{importGeneralModality}
#'   Runs a given formula on a subset of the data, and returns the model results. This is meant to help during the model selection process. \code{\link[glmmTMB]{glmmTMB}}. 
#'
#' @param ExperimentObj A SummarizedExperiment-type object generated from importGeneralModality
#' @param assayName a character string, matching the name of an assay within the SummarizedExperiment. 
#'   The assay named will be used for modeling. 
#' @param modelFormula The formula to use with glmmTMB, in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the ExperimentObj metadata.
#' @param ziFormula A formula for zero inflation. By default, zero-inflated modeling is turned off by setting ziFormula = ~ 0
#' @param zi_threshold A threshold above which the ziFormula is used. If zero-inflation is below this value (e.g. 0.1, or 1/10 values are zero), then those values are censored. This threshold is ignored if ziFormula = ~ 0
#' @param family distribution family parameter, passed to glmmTMB to describe the data's distribution.
#'     Default is normal (gaussian()). See  \link[glmmTMB]{glmmTMB}.
#' @param featureSample integer or vector of characters describing the features to test. 
#' If integer, stratified sampling will occur to identify features to test, with a minimum of 10. 
#' If a character vector, then all strings must align to features you want to test. 
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return modelList a list of outputs from glmmTMB
#'
#' @keywords pilot_modeling
#' @export
pilot_General <- function(ExperimentObj,
                      assayName = NULL,
                      modelFormula,
                      ziFormula = ~ 0,
                      zi_threshold = 0,
                      numCores = 1,
                      family = stats::gaussian(),
                      featureSample = 100,
                      verbose = FALSE) {

  if(isChAIObject(ExperimentObj, type = 'data', returnType = TRUE) != 'General'){

    stop('ExperimentObj is not a ChAI General Modality object, generated via importGeneralModality.')

  }

  if (length(assayName) > 1) {
    stop(
      "More than one assay was provided. ",
      "assayName must be length 1. To run over multiple assays/cell types, ",
      "combine the matrices into one and wrap them in a SummarizedExperiment. Then set assayName to ",
      "the name of that matrix in the obejct."
    )
  } else if (
    !assayName %in% names(SummarizedExperiment::assays(ExperimentObj))
  ) {
    stop("assayName was not found within ExperimentObj.")
  }

  names(SummarizedExperiment::assays(ExperimentObj))[names(SummarizedExperiment::assays(ExperimentObj)) == assayName] = 'counts'

  modelList <- .pilotModels_generic(SE_Object = ExperimentObj,
                        modelFormula = modelFormula,
                        ziFormula = ziFormula,
                        zi_threshold = zi_threshold,
                        family = family,
                        modality = 'General_Model',
                        verbose = verbose,
                        numCores = numCores,
                        featureSample = featureSample)
  return(modelList)

}



#' @title pilotModels_generic 
#'
#' @description \code{mpilotModels_generic} Internal generic funciton for modeling. 
#' @param SE_Object A SummarizedExperiment object generated from
#' @param modelFormula The formula for the continuous data that should be used within glmmTMB. It should be in the
#'   format (exp ~ factors). All factors must be found in column names of the SE_Object colData (i.e. sample metadata).
#'   See \link[glmmTMB]{glmmTMB}.
#' @param ziFormula The formula for the zero-inflated data that should be used within glmmTMB. It should be in the
#'   format ( ~ factors). All factors must be found in column names
#'   of SE_Object colData (i.e. sample metadata).
#' @param family distribution family parameter, passed to glmmTMB to describe the data's distribution.
#'     Default is normal (gaussian()). See  \link[glmmTMB]{glmmTMB}.
#' @param zi_threshold Zero-inflated threshold ( range = 0-1), representing the fraction of samples with zeros.
#'           At or above this threshold, the zero-inflated modeling kicks in.
#' @param  featureSample integer or vector of characters describing the features to test. 
#' If integer, stratified sampling will occur to identify features to test, with a minimum of 10. 
#' If a character vector, then all strings must align to features you want to test. 
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return a list of model objects from glmmTMB
#'
#' @noRd
#'
.pilotModels_generic <- function(SE_Object,
                        modelFormula = NULL,
                        ziFormula = NULL,
                        family = stats::gaussian(),
                        zi_threshold = 0,
                        modality = 'General_Model',
                        verbose = FALSE,
                        numCores = 1,
                        featureSample = 100) {

  Sample <- meanVal <- bin <- NULL

  modelingData <- SummarizedExperiment::assays(SE_Object)[[1]]
  MetaDF <- as.data.frame(SummarizedExperiment::colData(SE_Object))

  if(methods::is(modelFormula, 'character')){
      modelFormula= stats::as.formula(modelFormula)
  }else if(!methods::is(modelFormula, 'formula')){
      stop('Formula not Recognized')
  }
      
  if(methods::is(ziFormula, 'character')){
      ziFormula = stats::as.formula(ziFormula)
  }
      
  if(!methods::is(ziFormula, 'formula')){
      stop('Formula not Recognized')
  }

  if (!all(all.vars(modelFormula) %in% c("exp", colnames(MetaDF)))) {
    stop("Model formula is not in the correct format (exp ~ factors) or model factors are not found in sample-level metadata.")
  }

  if (!all(all.vars(ziFormula) %in% c(colnames(MetaDF))) & length(all.vars(ziFormula)) > 0) {
    stop("factors from the ziFormula were not found in the metadata.")
  }
  variableList <- c(all.vars(modelFormula)[all.vars(modelFormula) != "exp"], all.vars(ziFormula))
   
    if(!is.numeric(featureSample) & !all(featureSample %in% rownames(modelingData))){
      
        stop('featureSample was not recognized. It must either be the number of random samples to select, or a vector characters that align with features to test')
        
    }else if(all(featureSample %in% rownames(modelingData))){
     
        rowIndices <-featureSample
        
    }else if(is.numeric(featureSample)){
    
      samplingDF <- data.frame(meanVal = rowMeans(modelingData), rowName = rownames(modelingData))
      samplingDF <- dplyr::group_by(dplyr::mutate(samplingDF, bin = dplyr::ntile(meanVal, 10)), bin)
      rowIndices <- dplyr::sample_n(samplingDF, round(featureSample/10))
      rowIndices <- rowIndices$rowName
        
    }else{
    
        stop('featureSample was not recognized. This error should not be hit.')
        
    }

        
  MetaDF <- dplyr::filter(MetaDF, Sample %in% colnames(modelingData))

  ## Log transform the FragmentNumbers so as to stabilize the model. But only if FragNumber is in the model. Same for CellCounts.
  if(any(colnames(MetaDF) %in% c('FragNumber'))){
    MetaDF$rawFragNumber = MetaDF$FragNumber
    MetaDF$FragNumber = log10(MetaDF$FragNumber+1)
  }
  if(any(colnames(MetaDF) %in% c('CellCounts'))){
    MetaDF$rawCellCounts = MetaDF$CellCounts
    MetaDF$CellCounts = log10(MetaDF$CellCounts+1)
  }

  # Subset metadata to just the variables needed. This minimizes overhead for parallelization
  MetaDF <- MetaDF[, colnames(MetaDF) %in% c("Sample", variableList)]

  if (verbose) {
    message("Modeling results.")
  }

  if(modality == 'scATAC'){
    modelingData = log2(modelingData + 1)
  }
    
  # Make your clusters for efficient parallelization
  modelList <- pbapply::pblapply(X = rowIndices, function(x) {
  
    df <- data.frame(
      exp = as.numeric(modelingData[x, ]),
      MetaDF, stringsAsFactors = FALSE
    )
    tryCatch({
        
        ##
      if(modality != 'scATAC_Model'){
         modelRes <- glmmTMB::glmmTMB(exp ~ 1,
          ziformula = ~ 0,
          data = df,
          family = family,
          REML = TRUE,
          control = glmmTMB::glmmTMBControl(parallel = numCores)
        )
        
        sigZI <- DHARMa::testZeroInflation(modelRes)$p.value < 0.05
          
      }else{
      
          sigZI = FALSE
          
      }
        
     if(modality != 'scATAC_Model'){
         
        if(sigZI){
            modelRes <- glmmTMB::glmmTMB(modelFormula,
              ziformula = ziFormula,
              data = df,
              family = family,
              REML = TRUE,
              control = glmmTMB::glmmTMBControl(parallel = numCores)
             )
        }else{
          modelRes <- glmmTMB::glmmTMB(modelFormula,
              ziformula = ~ 0,
              data = df,
              family = family,
              REML = TRUE,
              control = glmmTMB::glmmTMBControl(parallel = numCores)
            )
         }
  
     }else if(sum(df$exp == 0, na.rm = TRUE) == 0){
        modelRes <- glmmTMB::glmmTMB(modelFormula,
          ziformula = ~ 0,
          data = df,
          family = family,
          REML = TRUE,
          control = glmmTMB::glmmTMBControl(parallel = numCores)
        )
      }else if(sum(df$exp == 0, na.rm = TRUE)/length(df$exp) <= zi_threshold){
        df$exp[df$exp == 0] = NA
        modelRes <- glmmTMB::glmmTMB(modelFormula,
          ziformula = ziFormula,
          data = df,
          family = family,
          REML = TRUE,
          control = glmmTMB::glmmTMBControl(parallel = numCores)
        )
        
        }else{
        modelRes <- glmmTMB::glmmTMB(modelFormula,
          ziformula = ziFormula,
          data = df,
          family = family,
          REML = TRUE,
          control = glmmTMB::glmmTMBControl(parallel = numCores)
        )
      }
    }, error = function(e){
      list('error' = e, 'Measurement' = x, 'Data' = df)
    })

  }, cl = NULL)
  names(modelList) <- rowIndices

  return(modelList)
}

