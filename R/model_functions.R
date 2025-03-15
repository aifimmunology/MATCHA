#' @title High-throughut Modeling of scATAC data
#'
#' @description \code{model_scATAC} High-throughput modeling of tile accessibility from scATAC \code{\link[glmmTMB]{glmmTMB}}. 
#'
#' @param atacSE A SummarizedExperiment object generated from
#'   getSampleTileMatrix. 
#' @param cellPopulation Name of a cell type. 
#' @param modelFormula The formula for the continuous data that should be used within glmmTMB. It should be in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the atacSE metadata, except for CellType, FragNumber and CellCount, which will be extracted from the atacSE.
#'   modelFormula must start with 'exp' as the response.
#'   See \link[glmmTMB]{glmmTMB}.
#' @param ziFormula The formula for the zero-inflated data that should be used within glmmTMB. It should be in the
#'   format ( ~ factors). All factors must be found in column names
#'   of the atacSE colData metadata, except for CellType, FragNumber and CellCount, which will be extracted from the atacSE.
#' @param zi_threshold Zero-inflated threshold ( range = 0-1), representing the fraction of samples with zeros. At or above this threshold, the zero-inflated modeling kicks in.
#' @param initialSampling Size of data to use for pilot
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return results a SummarizedExperiment containing ZI-GLMM results
#'
#'
#'
#' @examples
#' \dontrun{
#'   modelList <- model_scATAC(STM[c(1:1000),], 
#'                  cellPopulation = 'CD16 Mono',
#'                  modelFormula = exp~ Age + Sex + days_since_symptoms + (1|PTID), 
#'                  ziFormula = ~ 0 + FragNumber + Age, 
#'                  verbose = TRUE, 
#'                  numCores = 35 )
#' }
#'
#' @export
#' @keywords modeling_individual
model_scATAC <- function(atacSE,
                      cellPopulation,
                      modelFormula = NULL,
                      ziFormula = ~ 0,
                      zi_threshold = 0,
                      initialSampling = 5,
                      verbose = FALSE,
                      numCores = 2) {
  if (!requireNamespace("MOCHA", quietly = TRUE)) {
      stop(
      "Package 'MOCHA' is required for pilot_scATAC. ",
      "Please install 'MOCHA' to proceed."
      )
  }

  if(isChAIObject(atacSE, type = 'data', returnType = TRUE) != 'scATAC'){

    stop('atacSE is not a MOCHA SampleTileObject.')

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
      "use MOCHA to run combineSampleTileMatrix() to produce a new combined atacSE and set ",
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

  exp <- .model_generic(SE_Object = newObj,
                      continuousFormula = modelFormula,
                      ziFormula = ziFormula,
                      zi_threshold = zi_threshold,
                      initialSampling = initialSampling,
                      family = stats::gaussian(),
                      modality = 'scATAC_Model',
                      verbose = verbose,
                      numCores = numCores)
  return(exp)

}


#' @title High-throughout Modeling of Gene Expression 
#'
#' @description \code{model_scRNA} Runs high-throughout GLM modeling of gene expression using \code{\link[glmmTMB]{glmmTMB}}
#'
#' @param rnaSE A SummarizedExperiment object generated from normalizePseudobulk
#' @param cellPopulation Name of a cell type. 
#' @param modelFormula The formula to use, in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the ExperimentObj metadata. modelFormula must start with 'exp' as the response.
#'   See \link[glmmTMB]{glmmTMB}.
#' @param ziFormula A formula for zero inflation. By default, zero-inflated modeling is turned off by setting ziFormula = ~ 0
#' @param family String. Can be 'negativeBinomial1', 'negativeBinomial2', or 'poisson'. Default is "negativeBinomial2". 
#' @param detectionThreshold A number between 0 and 1, representing the mean detection rate threshold for a given gene to be modeled. 
#'   This detection rate is calculated for each sample and cell type during makePseudobulkRNA, and represents the percentage of cells that have a transcript for a given gene. Over all samples, the average has to be above this to be modeled. Default is 0.01. 
#' @param expressionThreshold A number greater than zero, representing the expression threshold for modeling. A given gene, on average across all samples, expressed above this threshold. The default is 0. 
#' @param cellCountThreshold The minimum number of cells in a given pseudobulk for the pseudobulk to be included in analysis. If fewer than this number of cells are found, then the sample will be dicarded The number of cells within the pseudobulked scRNA. Default is 10 cells. 
#' @param initialSampling Size of data to use for pilot
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return results a SummarizedExperiment containing LMEM results. Assays are metrics related to the model coefficients,
#'          including the Estimate, Std_Error, df, t_value, p_value. Within each assay, each row corresponds to each row of
#'          the SummarizedExperiment and columns correspond to each fixed effect variable within the model.
#'          Any row metadata from the ExperimentObject (see rowData(ExperimentObj)) is preserved in the output. 
#'          The Residual matrix and the variance of the random effects are saved in the metadata slot of the output. 
#'
#'
#'
#' @examples
#' \dontrun{
#'   modelList <- model_scRNA(rnaSE,
#'    cellPopulation = "CD14_Mono"
#'     modelFormula = "exp ~ treatment + (1|patientID)",
#'     initialSampling = 5,
#'     verbose = FALSE,
#'     numCores = 1
#'  )
#' }
#'
#' @export
#' @keywords modeling_individual
model_scRNA <- function(rnaSE,
                    cellPopulation = NULL,
                    modelFormula = NULL,
                    ziFormula = ~0,
                    zi_threshold = 0,
                    family = "negativeBinomial2",
                    detectionThreshold = 0.01,
                    expressionThreshold = 0,
                    cellCountThreshold = 10,
                    initialSampling = 5,
                    verbose = FALSE,
                    numCores = 2) {

  if(isChAIObject(rnaSE, type = 'data', returnType = TRUE) != 'scRNA'){

    stop('rnaSE is not a ChAI scRNA Object (normalized, pseudobulked via ChAI.)')

  }
  
  if(!any(cellPopulation %in% names(SummarizedExperiment::assays(rnaSE)))){
    stop('ExperimentObj does not contain an assay that matches the assayName input variable.')
  }
  if(length(cellPopulation) != 1){
    stop('Please provide only one cell type at time, or run flattenChAI to run over multiple celltypes.')
  }
  # Subset object to just one cell type, and then combine to flatten additional metadata (like cell counts) into the summarizedExperiment metadata. 
  newRNA = flattenChAI(rnaSE,  cellPopulations = cellPopulation, metadataT = TRUE)
  
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
      
  if (!methods::is(modelFormula,'formula')) {
    stop("modelFormula was not provided as a formula.")
  }
    
    if (!methods::is(ziFormula,'formula')) {
    stop("ziFormula was not provided as a formula.")
  }

  if (zi_threshold < 0 | zi_threshold > 1 | ! is.numeric(zi_threshold)) {
    stop("zi_threshold must be between 0 and 1.")
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

          message(stringr::str_interp('${sum(!newRNA$CellCounts >=  cellCountThreshold)} Samples removed due to cell counts below cellCountThreshold'))

          newRNA <- newRNA[,newRNA$CellCounts >=  cellCountThreshold]

      }
     
  if(verbose){
    message(stringr::str_interp("${NROW(newRNA)} genes will be modeled")) 
  }
    
 
  
  exp <- .model_generic(SE_Object = newRNA,
                      continuousFormula = stats::as.formula(modelFormula),
                      ziFormula = ziFormula,
                      zi_threshold = zi_threshold,
                      family = family,
                      initialSampling = initialSampling,
                      modality = 'scRNA_Model',
                      verbose = verbose,
                      numCores = numCores)
  return(exp)

}


#' @title High-through modeling (Generalized)
#'
#' @description \code{model_General} Runs generalized GLM modeling for
#'   continuous, non-zero inflated data using \code{\link[glmmTMB]{glmmTMB}}
#'
#' @param ExperimentObj A SummarizedExperiment object generated from chromVAR, or other.
#'   It is expected to contain normally distributed data without zero-inflation. 
#' @param assayName The name of the assay to model within the SummarizedExperiment. 
#' @param modelFormula The formula to use, in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the ExperimentObj metadata. modelFormula must start with 'exp' as the response.
#'   See \link[glmmTMB]{glmmTMB}.
#' @param ziFormula A formula for zero inflation. By default, zero-inflated modeling is turned off by setting ziFormula = ~ 0
#' @param family distribution family parameter, passed to glmmTMB to describe the data's distribution.
#'     Default is normal (gaussian()). See  \link[glmmTMB]{glmmTMB}.
#' @param initialSampling Size of data to use for pilot
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return results a SummarizedExperiment containing LMEM results. Assays are metrics related to the model coefficients,
#'          including the Estimate, Std_Error, df, t_value, p_value. Within each assay, each row corresponds to each row of
#'          the SummarizedExperiment and columns correspond to each fixed effect variable within the model.
#'          Any row metadata from the ExperimentObject (see rowData(ExperimentObj)) is preserved in the output. 
#'          The Residual matrix and the variance of the random effects are saved in the metadata slot of the output. 
#'
#'
#'
#' @examples
#' \dontrun{
#'   modelList <- model_General(ExperimentObj,
#'    assayName = "z",
#'     modelFormula = NULL,
#'     initialSampling = 5,
#'     verbose = FALSE,
#'     numCores = 1
#'  )
#' }
#'
#' @export
#' @keywords modeling_individual
model_General <- function(ExperimentObj,
                    assayName,
                    modelFormula,
                    ziFormula = ~ 0,
                    family = stats::gaussian(),
                    initialSampling = 5,
                    verbose = FALSE,
                    numCores = 2) {

  if(isChAIObject(ExperimentObj, type = 'data', returnType = TRUE) != 'General'){

    stop('ExperimentObj is not a ChAI General Modality object, generated via importGeneralModality.')

  }
    
   if (methods::is(ziFormula, 'character')) {
    ziFormula <- stats::as.formula(ziFormula)
  }

  
  if(!any(names(SummarizedExperiment::assays(ExperimentObj)) %in% assayName)){
    stop('ExperimentObj does not contain an assay that matches the assayName input variable.')
  }
  SummarizedExperiment::assays(ExperimentObj) = SummarizedExperiment::assays(ExperimentObj)[assayName]

  exp <- .model_generic(SE_Object = ExperimentObj,
                      continuousFormula = stats::as.formula(modelFormula),
                      ziFormula = ziFormula,
                      zi_threshold = 0,
                      family = family,
                      initialSampling = initialSampling,
                      modality = 'General_Model',
                      verbose = verbose,
                      numCores = numCores)
  return(exp)

}




#' @title High-throughput Modeling of ChromVAR Motif scores
#'
#' @description \code{model_ChromVAR} Runs linear mixed-effects modeling for
#'   chromVAR Z-scores or deviations, derived from a MOCHA object using makeChromVAR, and modeled using \code{\link[glmmTMB]{glmmTMB}}
#'
#' @param chromSE A SummarizedExperiment object generated from chromVAR via ChAI::makeChromVAR
#' @param cellPopulation The name of the assay to model within the SummarizedExperiment. 
#' @param modelFormula The formula to use, in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the chromSE metadata. modelFormula must start with 'exp' as the response.
#'   See \link[glmmTMB]{glmmTMB}.
#' @param initialSampling Size of data to use for pilot
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return results a SummarizedExperiment containing LMEM results. Assays are metrics related to the model coefficients,
#'          including the Estimate, Std_Error, df, t_value, p_value. Within each assay, each row corresponds to each row of
#'          the SummarizedExperiment and columns correspond to each fixed effect variable within the model.
#'          Any row metadata from the chromSE (see rowData(chromSE)) is preserved in the output. 
#'          The Residual matrix and the variance of the random effects are saved in the metadata slot of the output. 
#'
#'
#'
#' @examples
#' \dontrun{
#'   modelList <- model_ChromVAR(chromSE,
#'    cellPopulation = "CD16 Mono",
#'     modelFormula = NULL,
#'     initialSampling = 5,
#'     verbose = FALSE,
#'     numCores = 1
#'  )
#' }
#'
#' @export
#' @keywords modeling_individual
#' 

model_ChromVAR <- function(chromSE,
                    cellPopulation,
                    modelFormula,
                    initialSampling = 5,
                    verbose = FALSE,
                    numCores = 2) {

  if(isChAIObject(chromSE, type = 'data', returnType = TRUE) != 'ChromVAR'){

    stop('chromSE is not a ChAI ChromVAR Object, created via makeChromVAR.')

  }
  
  if(!any(names(SummarizedExperiment::assays(chromSE)) %in% cellPopulation)){
    stop('chromSE does not contain an assay that matches the cellPopulation input variable.')
  }
  if(length(cellPopulation) != 1){
    stop('cellPopulation must be only one string, matching a cell population within the ChAI-chromVAR object.')
  }
  # Subset object to just one cell type, and then combine to flatten additional metadata (like cell counts) into the summarizedExperiment metadata. 
  newChrom = flattenChAI(chromSE, cellPopulations = cellPopulation)

  exp <- .model_generic(SE_Object = newChrom,
                      continuousFormula = stats::as.formula(modelFormula),
                      ziFormula = ~ 0,
                      zi_threshold = 0,
                      family = stats::gaussian(),
                      initialSampling = initialSampling,
                      modality = 'ChromVAR_Model',
                      verbose = verbose,
                      numCores = numCores)
  return(exp)

}





#' @title model_generic 
#'
#' @description \code{model_generic} Internal generic funciton for modeling. 
#' @param SE_Object A SummarizedExperiment object generated from
#' @param continuousFormula The formula for the continuous data that should be used within glmmTMB. It should be in the
#'   format (exp ~ factors). All factors must be found in column names of the SE_Object colData (i.e. sample metadata).
#'   See \link[glmmTMB]{glmmTMB}.
#' @param ziFormula The formula for the zero-inflated data that should be used within glmmTMB. It should be in the
#'   format ( ~ factors). All factors must be found in column names
#'   of SE_Object colData (i.e. sample metadata).
#' @param zi_threshold Zero-inflated threshold ( range = 0-1), representing the fraction of samples with zeros.
#'           At or above this threshold, the zero-inflated modeling kicks in.
#' @param initialSampling Size of data to use for pilot
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return results a SummarizedExperiment containing ZI-GLMM results
#'
#' @noRd
#'
.model_generic <- function(SE_Object,
                      continuousFormula = NULL,
                      ziFormula = NULL,
                      zi_threshold = 0,
                      initialSampling = 5,
                      family = stats::gaussian(),
                      modality = 'General_Model',
                      verbose = FALSE,
                      numCores = 2) {
  Sample <- NULL
                 
  if (!methods::is(continuousFormula, "formula")) {

    stop("continuousFormula was not provided as a formula.")

  }
  if (!methods::is(ziFormula, "formula") & !is.null(ziFormula)) {

    stop("ziFormula was not provided as a formula.")

  }

  if(modality == 'scATAC_Model'){

    modelingData <- log2(SummarizedExperiment::assays(SE_Object)[[1]]+1)

  }else if(modality == 'scRNA_Model'){

    modelingData <- SummarizedExperiment::assays(SE_Object)[[1]]

  }else{

    modelingData <- SummarizedExperiment::assays(SE_Object)[[1]]

  }

  modelingData[is.na(modelingData)]= 0 ## We need this for ChromVAR where a sample can fall if too few reads are present. In future, we should maybe make this edits specific to chromVAR and through an error for other modalities. scATAC replaces NAs with zeros during combineSampleTileMatrix
    
  #Remove any samples that had no cells detected 
  emptySamples <- unlist(apply(modelingData, 2, function(x) all(x == 0)))
                               
  MetaDF <- as.data.frame(SummarizedExperiment::colData(SE_Object))
                    
  if(any(emptySamples)){
    tmpNum <- sum(emptySamples)
    warning(stringr::str_interp("${tmpNum} empty/missing samples found (only zero values). Removing these samples. This can occur in scATAC/scRNA when a given sample has no cells of a given cell type."))
    modelingData[,emptySamples] = NA
  }

  if (!all(all.vars(continuousFormula) %in% c("exp", colnames(MetaDF)))) {
    stop("Model formula is not in the correct format (exp ~ factors) or model factors are not found in column names of metadata.")
  }

  if (!all(all.vars(ziFormula) %in% colnames(MetaDF))) {
    stop("factors from the ziFormula were not found in the metadata.")
  }

  if (!"Sample" %in% colnames(MetaDF)) {
    stop("There is no Sample column within the metadata of the object you are modeling. Please add a Sample column to the SummarizedExperiment obj and try again.")
  }

  if (!all(MetaDF$Sample  == colnames(modelingData))) {
    stop("The Sample metadata column of the SummarizedExperiment do not align with the SummarizedExperiment's sample names. Please check the Sample column, or generate a new SummarizedExperiment object using ChAI::importGeneralModality.")
  }

  variableList <- c(all.vars(continuousFormula)[all.vars(continuousFormula) != "exp"], all.vars(ziFormula))

  # Subset metadata to just the variables needed. This minimizes overhead for parallelization.
  # The Sample column is needed for alignment between data.frames later. 
  MetaDF <- MetaDF[, colnames(MetaDF) %in% c('Sample',variableList)]

  #Transform in character strings for multithreading. 
  continuousFormula <- deparse(continuousFormula)
  ziFormula <- deparse(ziFormula)

  ## Log transform the FragmentNumbers so as to stabilize the model. But only if FragNumber is in the model. Same for CellCounts.
  if(any(colnames(MetaDF) %in% c('FragNumber'))){
    MetaDF$rawFragNumber = MetaDF$FragNumber
    MetaDF$FragNumber <- log10(MetaDF$FragNumber+1)
  }
  if(any(colnames(MetaDF) %in% c('CellCounts'))){
    MetaDF$rawCellCounts = MetaDF$CellCounts
    MetaDF$CellCounts <- log10(MetaDF$CellCounts+1)
  }

  if (verbose) {
    message("Running a quick test.")
  }


  # Generate pilot data for the null data.frame
  nullDFList <- generateNULL(modelingData, MetaDF,continuousFormula, ziFormula, family, modality, initialSampling)
  if (verbose) {
    message("Modeling results.")
  }

  if(numCores <= 1){
    stop('numCores must be greater than 1. This method is meant to be parallelized.')
  }
  # Make your clusters for efficient parallelization

  cl <- parallel::makeCluster(numCores)

  iterList = lapply(rownames(modelingData), function(x){

           data.frame(
                  exp = as.numeric(modelingData[x, ]),
                  MetaDF, stringsAsFactors = FALSE
                )
  })
                           
  suppressMessages(parallel::clusterEvalQ(cl, {
        suppressMessages(library('glmmTMB'))
      suppressMessages(library('DHARMa'))
      }))
    

                               
  coeffList <- pbapply::pblapply(cl = cl, X = iterList, individualZIGLMM, 
                                continuousFormula1 = continuousFormula,
                                ziFormula1 = ziFormula,
                                family1 = family,
                                zi_threshold1 = zi_threshold,
                                nullDFList1 = nullDFList,
                                modality1 = modality)
                               
  parallel::stopCluster(cl)
  if (verbose) {
      message("Reorganizing residuals and random effect variance.")
  }

  if(!all(is.null(SummarizedExperiment::rowRanges(SE_Object)))){
    ranged = TRUE
  }else{
    ranged= FALSE
  }

  processedOuts <- processModelOutputs(modelOutputList = coeffList, 
                                        nullDFList = nullDFList, 
                                        rownamesList = rownames(modelingData),
                                        ranged = ranged,
                                        SummarizedExperimentObj = SE_Object
                                        )
  processedOuts@metadata = append(processedOuts@metadata, list('Type' = modality))
  return(processedOuts)
}

#' @title Internal function to run generating null results, in case of model failure
#'
#' @description \code{generateNULL} Runs null results, either for glmmTMB
#' @param modelingData the data used for modeling. 
#' @param MetaDF the metadata associated with the data
#' @param continuousFormula The formula for the continuous data that should be used within glmmTMB. It should be in the
#'   format (exp ~ factors). All factors must be found in column names of the SE_Object colData (i.e. sample metadata).
#'   See \link[glmmTMB]{glmmTMB}.
#' @param ziFormula The formula for the zero-inflated data that should be used within glmmTMB. It should be in the
#'   format ( ~ factors). All factors must be found in column names of MetaDF
#' @param family family for data distribution to be used in the model
#' @param modality a flag to mark whether the data is from scATAC, scRNA, or General. 
#' @param initialSampling Size of data to use for pilot
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @return output_vector A linear model
#'
#' @noRd

generateNULL <- function(modelingData, MetaDF, continuousFormula, ziFormula, family, modality, initialSampling){

  pilotIndices <- sample(x = 1:dim(modelingData)[1], size = initialSampling, replace = FALSE)
  message('Sampling data. Identifying null template.')
  modelList <- pbapply::pblapply(X = pilotIndices, function(x) {
    df <- data.frame(
      exp = as.numeric(modelingData[x, ]),
      MetaDF, stringsAsFactors = FALSE
    )

    tryCatch(
      {
        suppressWarnings(glmmTMB::glmmTMB(stats::as.formula(paste(continuousFormula, collapse = " ")),
            ziformula = stats::as.formula(paste(ziFormula, collapse = " ")),
            data = df,
            family = family,
            REML = TRUE))
       
      },
      error = function(e) {
        NA
      }
    )
  }, cl = NULL)

  NAList <- unlist(lapply(modelList, function(x){
                    ## First identify all the pilot models that ran at all without errors. 
                   tmp1 <- all(is.na(x))
                   if(tmp1){ 
                    return(FALSE)
                   }else{
                     ## Then identify all the models that converged. 
                    tryCatch(
                      {
                        tmp1 <- summary(x)$coefficients
                        return(TRUE)
                      }, error = function(e){ return(FALSE)})
                   }
            }))
   
  if (all(!NAList)) {
    stop("For the initial sampling, every test model failed. Reconsider modelFormula or increase initial sampling size.")
  } else {
    idx <- which(NAList)
    # Did any of the models have both zero-inflated and continous portions?
    bothZI_Cont <- unlist(lapply(idx, function(x){
      ziDF <- summary(modelList[[x]])$coefficients$zi
      !is.null(ziDF) | sum(dim(ziDF)) != 0

    }))

    if(all(!bothZI_Cont) & !ziFormula %in% c('~ 0', '~0')){
      warning('No working models using the zero-inflated formula. Do you need to modify the zero-inflated formula?')
      modelRes <- modelList[[idx[1]]]
    }else if(ziFormula %in% c('~ 0', '~0')){
      modelRes <- modelList[[idx[1]]]
    }else {
      modelRes <- modelList[[idx[which(bothZI_Cont)[1]]]]
    }

    #Extract the first representative model. And use it to create a null template. 
    coeff2 <- summary(modelRes)$coefficients
    coeff <- lapply(coeff2, as.data.frame)
    coeff$cond[!is.na(coeff$cond)] <- NA
    rownames(coeff$cond)[grepl('(Intercept)',rownames(coeff$cond))] = 'Intercept'
    if(all(!bothZI_Cont)){
      combinedCoeff <- coeff$cond
    }else{
      coeff$zi[!is.na(coeff$zi)] <- NA
      rownames(coeff$zi)[grepl('(Intercept)',rownames(coeff$zi))] = 'Intercept'
      rownames(coeff$zi) <- paste('ZI', rownames(coeff$zi), sep ='_')
      combinedCoeff <- rbind(coeff$cond, coeff$zi)
    }
    
    Resid <- stats::resid(modelRes)
    
    if(!all(MetaDF$Sample %in% names(Resid))){
      NA_samples <- rep(NA, sum(!MetaDF$Sample %in% names(Resid)))
      names(NA_samples) <- MetaDF$Sample[!MetaDF$Sample %in% names(Resid)]
      Resid <- c(Resid, NA_samples)
    }
    Resid <- Resid[match(names(Resid),MetaDF$Sample)]
    Resid[!is.na(Resid)] <- NA
    varCorrObj <- glmmTMB::VarCorr(modelRes)
    if(is.null(varCorrObj$cond) & all(!bothZI_Cont)){
      varCorrObj <- list('cond' = list('NoRandomEffects'= NA), 'zi' = NULL)
      residual = c('Residual' = NA)
    } else if(is.null(varCorrObj$cond)){
      varCorrObj <- list('cond' = list('NoRandomEffects'= NA), 'zi' = list('None'= NA))
      residual = c('Residual' = NA)
    }else{
      residual = as.vector(attr(varCorrObj$cond, "sc")^2)
      names(residual) = 'Residual'
    }
    cond_other = unlist(varCorrObj$cond)
    names(cond_other) = paste('Cond', names(cond_other), sep = "_")
    
   
    if(!is.null(varCorrObj$zi)){
      zi_other = unlist(varCorrObj$zi)
      names(zi_other) = paste('ZI', names(zi_other), sep = "_")
      varcor_df <- c(cond_other, zi_other,residual)
    }else{
      varcor_df <- c(cond_other, residual)
    }
    varcor_df[!is.na(varcor_df)] = NA

    nullDFList <- list('Coeff' = combinedCoeff, 'Resid' = Resid, 'VCov'= varcor_df)

    rm(modelList)
  }
  return(nullDFList)
}

#' @title \code{IndividualZIGLMM}
#'
#' @description \code{IndividualZIGLMM} Runs zero-inflated linear modeling on data provided. Written for efficient parallelization.
#'
#' @param refList. A list where the first index is a data.frame to use for modeling, and the second is the formula for modeling.
#'
#' @return A linear model
#'
#' @noRd
#'


individualZIGLMM <- function(df, continuousFormula1, family1, ziFormula1, zi_threshold1, nullDFList1, modality1) {

    ## Test for zero-inflation, if it isn't a scATAC modeling question. 
   sigZI <- tryCatch(
    {
     if(modality1 != 'scATAC_Model' & all(ziFormula != '~ 0')){
         suppressWarnings(modelRes2 <- glmmTMB::glmmTMB(formula= exp ~ 1,
          ziformula = ~ 0,
          data = df,
          family = family1,
          REML = TRUE,
          control = glmmTMB::glmmTMBControl(parallel = 1)
        ))
        
        tmpZI <- suppressWarnings(DHARMa::testZeroInflation(modelRes2)$p.value < 0.05)
        if(is.na(tmpZI)){  
            FALSE
        }else{
            tmpZI
        }
         
          
      }else{
      
        FALSE
          
      }
    },
    error = function(e) {
       FALSE
   }
   )
        
   output_vector <- tryCatch(
    {
      # for modeling that isn't scATAC, run zero-inflated modeling only on zero-inflated features. 
     if(modality1 != 'scATAC_Model'){
         
        if(sigZI){
            suppressWarnings(modelRes <- glmmTMB::glmmTMB(stats::as.formula(paste(continuousFormula1, collapse = " ")),
                ziformula = stats::as.formula(paste(ziFormula1, collapse = " ")),
              data = df,
              family = family1,
              REML = TRUE,
              control = glmmTMB::glmmTMBControl(parallel = 1)
             ))
            
        }else{
         suppressWarnings( modelRes <- glmmTMB::glmmTMB(stats::as.formula(paste(continuousFormula1, collapse = " ")),
              ziformula = ~ 0,
              data = df,
              family = family1,
              REML = TRUE,
              control = glmmTMB::glmmTMBControl(parallel = 1)
            ))
         }
  
     }else if(sum(df$exp == 0, na.rm = TRUE)/length(df$exp) == 0){
        suppressWarnings(modelRes <- glmmTMB::glmmTMB(stats::as.formula(paste(continuousFormula1, collapse = " ")),
          ziformula = ~ 0,
          data = df,
          family = family1,
          REML = TRUE,
          control = glmmTMB::glmmTMBControl(parallel = 1)
        ))
       
      }else if(sum(df$exp == 0, na.rm = TRUE)/length(df$exp) <= zi_threshold1){
        df$exp[df$exp == 0] = NA

        suppressWarnings(modelRes <- glmmTMB::glmmTMB(stats::as.formula(paste(continuousFormula1, collapse = " ")),
          ziformula = ~ 0,
          data = df,
          family = family1,
          REML = TRUE,
          control = glmmTMB::glmmTMBControl(parallel = 1)
        ))
       
      }else {

        suppressWarnings(modelRes <- glmmTMB::glmmTMB(stats::as.formula(paste(continuousFormula1, collapse = " ")),
           ziformula = stats::as.formula(paste(ziFormula1, collapse = " ")),
          data = df,
          family = family1,
          REML = TRUE,
          control = glmmTMB::glmmTMBControl(parallel = 1)
        ))
      }

      Coeff <- lapply(summary(modelRes)$coefficients, as.data.frame)
      rownames(Coeff[['cond']])[grepl('(Intercept)',rownames(Coeff[['cond']]))] = 'Intercept'
      if(sum(dim(Coeff$zi)) != 0){
        rownames(Coeff$zi)[grepl('(Intercept)',rownames(Coeff$zi))] = 'Intercept'
        rownames(Coeff$zi) <- paste('ZI', rownames(Coeff$zi), sep ='_')
      }else{
        Coeff$zi = nullDFList1$Coeff[grepl('ZI',rownames(nullDFList1$Coeff)),]
      }
      combinedCoeff <- rbind(Coeff$cond, Coeff$zi)


      Resid <- stats::resid(modelRes)
      
      #Clean up residuals. Data that is NA will be removed from the residuals, so we add them back in as NAs for the sake of completion. 
      if(!all(df$Sample %in% names(Resid))){
        NA_samples <- rep(NA, sum(!df$Sample %in% names(Resid)))
        names(NA_samples) <- df$Sample[!df$Sample %in% names(Resid)]
        Resid <- c(Resid, NA_samples)
      }
      Resid <- Resid[match(names(Resid),df$Sample)]

      #Now process the variance from random effects. 
      varcor_df <- tryCatch({
        varCorrObj <- glmmTMB::VarCorr(modelRes)
        cond_other = unlist(varCorrObj$cond)
        names(cond_other) = paste('Cond', names(cond_other), sep = "_")
        residual = as.vector(attr(varCorrObj$cond, "sc")^2)
        names(residual) = 'Residual'

        #Process variance
        if(!is.null(varCorrObj$zi)){
          zi_other = unlist(varCorrObj$zi)
          names(zi_other) = paste('ZI', names(zi_other), sep = "_")
          varcor_df <- c(cond_other, zi_other,residual)
        }else if(all(df$exp !=0, na.rm = TRUE)){
          subNull = nullDFList[grepl('ZI_', names(nullDFList))]
          zi_other = rep(0, length(subNull))
          names(zi_other) = names(subNull)
          varcor_df <- c(cond_other, zi_other,residual)
        }else {
          varcor_df <- c(cond_other, residual)
        }
        varcor_df
      }, error = function(e){
        nullDFList1$VCov
      })
      
      return(list('Coeff' = combinedCoeff, 'Resid' = Resid, 'VCov'= varcor_df))
    },
    error = function(e) {
      nullDFList1
    }
  )
  return(output_vector)
}



#' @title Internal function to processing model outputs
#'
#' @description \code{processModelOutputs} 
#' @param modelOutputList. A list of modeloutputs, processed by model_scATAC, model_scRNA, or model_General. 
#'        The first output is the coefficient data.frame, then the residuals, and the Variance. 
#' @param nullDFList A null templates for the model outputs
#' @param rownamesList a name of all the rows that were interate over. 
#' @param ranged A boolean, with default FALSE. This determines whether the data is involves GenomicRanges (like RNA - gene bodies, or ATAC- tiles).
#' @param SummarizedExperimentObj SummarizedExperiment object used for modeling. From this, rowData and colData will be preserved with the residuals.
#' @param returnList A boolean, default is FALSE. For associative models, returnList = TRUE in order to return the values without repackaging into a SummarizedExperiment.
#' @return A SummarizedExperiment object, that captures the models performance. Each fixed effect will be one assay, columns will be the 
#'        statistics for the fixed effect and measurements (Estimate, Error, p-value, etc..). Residuals and Variance will be saved in the object's metadata.
#'
#' @noRd
processModelOutputs <- function(modelOutputList, nullDFList, rownamesList, ranged = FALSE,
                                  SummarizedExperimentObj, returnList = FALSE) {
   
    coeffNames <- rownames(nullDFList$Coeff)
    newColumnNames <- gsub('Pr\\(>\\|.\\|)','p_value', gsub(' |\\. ','_',colnames(nullDFList$Coeff)))
    output_list <- lapply(coeffNames, function(z){
      tmpCoef <- do.call("rbind", lapply(X = modelOutputList, function(x) {
            tmpDf <- x[['Coeff']][z,]
            colnames(tmpDf) <- newColumnNames
            tmpDf
          }))
      rownames(tmpCoef) <- rownamesList
      tmpCoef$FDR <- stats::p.adjust(tmpCoef$p_value, 'fdr')
      tmpCoef
    })
    names(output_list) <- gsub('Pr\\(>\\|.\\|)','p_value', gsub(' |\\. ','_',coeffNames))
    if(grepl('Intercept', names(output_list)[1])){
      names(output_list)[1] = 'Intercept'
    }

    residual_tmp <- do.call(
      "rbind", lapply(X = modelOutputList, function(x) {
        x[['Resid']]
      })
    )
    vcov_tmp <- do.call(
      "rbind", lapply(X = modelOutputList, function(x) {
        x[['VCov']]
      })
    )
    rownames(residual_tmp) <- rownames(vcov_tmp) <- rownamesList
    residual_tmp <- residual_tmp[,match(rownames(SummarizedExperiment::colData(SummarizedExperimentObj)), colnames(residual_tmp))]

    if(returnList){
      return(list('output' = output_list , 'Resid' = residual_tmp , 'Variance' = vcov_tmp))
    }
    
    #Repackage Residuals into a SummarizedExperiment
    if(ranged){
      ResidualSE <- SummarizedExperiment::SummarizedExperiment(
                      list('Residual' =  residual_tmp),
                      colData = SummarizedExperiment::colData(SummarizedExperimentObj),
                      rowRanges = SummarizedExperiment::rowRanges(SummarizedExperimentObj),
                      metadata = S4Vectors::metadata(SummarizedExperimentObj)
              )
      #Package up the metadata list. 
      metaDataList <- list('Residuals' = ResidualSE,
                  'RandomEffectVariance' = vcov_tmp)

      results <- SummarizedExperiment::SummarizedExperiment(
          output_list,
          rowRanges = SummarizedExperiment::rowRanges(SummarizedExperimentObj),
          metadata = metaDataList
      )
    }else{
      ResidualSE <- SummarizedExperiment::SummarizedExperiment(
                  list('Residual' =  residual_tmp),
                  colData = SummarizedExperiment::colData(SummarizedExperimentObj),
                  rowData = SummarizedExperiment::rowData(SummarizedExperimentObj),
                  metadata = S4Vectors::metadata(SummarizedExperimentObj)
          )
      #Package up the metadata list. 
      metaDataList <- list('Residuals' = ResidualSE,
                        'RandomEffectVariance' = vcov_tmp)
     results <- SummarizedExperiment::SummarizedExperiment(
       output_list,
       rowData = SummarizedExperiment::rowData(SummarizedExperimentObj),
       metadata = metaDataList
     )
    
    }

    return(results)
}

