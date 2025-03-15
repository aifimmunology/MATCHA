#' @title \code{visualizeAssociations}
#'
#' @description \code{visualizeAssociations} Takes the output of multiModalModeling (or other) and visualizes the top interaction features
#'  for one of the modalities. 
#' @param interAct A dataframe of coefficients related to a modality, extracted from multiModalModeling or other. This should be prefiltered for significance (FDR < 0.1)
#' @param dataType1 Boolean flag. This determines whether it should look at the Modality1 or Modality2 column (i.e. the first or second assay given to multiModalIntegration)
#' @param dataTypeName A string to label the modality in the plots generated. 
#' @param threshold A number, describing how many points to label. By default, it will label the top 20 points. 
#' @param max.overlaps The maximum number of overlapping labels, a parameter which is passed to ggrepel when labeling points on the plot. 
#' @param verbose Boolean flag to determine verbosity. 
#' @param rastr A boolean flag to determine whether to rasterize the pdfs. Requires ggrastr to be installed.
#'
#' @return A ggplot object, showing number of interactions on the y-axis and the rank order of measurements on the x-axis 
#'    (Measurements rank-ordered by the number of interactions they have with the other modality.)
#'
#'  
#'
#'
#' @export
#' @keywords model_results

visualizeAssociations <- function(interAct, dataType1 = TRUE, dataTypeName = 'TF', threshold = 20,  max.overlaps = 10,
                                 verbose = FALSE, rastr = FALSE){
    Interactions <- Factor <- Rank <- NULL
    if(dataType1){
        factor1 = 'Modality1'
    }else{
        factor1 = 'Modality2'
    }
    df1 <- as.data.frame(table(interAct[,factor1]))
    colnames(df1) = c(dataTypeName,'Interactions')
    df1 <- df1[order(df1$Interactions,decreasing =T),]
    df1$Rank <- as.numeric(c(1:length(df1$Interactions)))
    df1$Factor <- df1[,dataTypeName]
    df1$Factor[df1$Rank > threshold] = NA
    if(verbose){
       print(df1[c(1:5),])
    }

    if(rastr){
        if (!requireNamespace("ggrastr", quietly = TRUE)) {
            stop(
                "Package 'ggrastr' is required for visualizeAssociations when rastr = TRUE. ",
                "Please install 'ggrastr' or set rastr = FALSE to proceed."
            )
        }
        ggplot2::ggplot(df1, ggplot2::aes(x = Rank, y = Interactions, label = Factor)) + 
        ggrastr::geom_point_rast() + ggplot2::theme_minimal() +
        ggplot2::ggtitle(factor1) + 
        ggplot2::xlab('Rank') + ggplot2::ylab('Number of Interactions')  + 
        ggrepel::geom_text_repel(max.overlaps = max.overlaps)
    }else{
        ggplot2::ggplot(df1, ggplot2::aes(x = Rank, y = Interactions, label = Factor)) + 
        ggplot2::geom_point() + ggplot2::theme_minimal() +
        ggplot2::ggtitle(factor1) + 
        ggplot2::xlab('Rank') + ggplot2::ylab('Number of Interactions')  + 
        ggrepel::geom_text_repel(max.overlaps = max.overlaps)
    }
}

#' @title \code{volcanoPairs}
#'
#' @description \code{volcanoPairs} Takes the output of multiModalModeling (or other) and generates a volcano plot of Estimate vs -log10(FDR) for that pair
#' @param interModel The summarizedExperiment output of any association modeling funciton (scATAC_Associations, scRNA_Associations, general_Associations, & GeneTile_Associations)
#' @param FDR_threshold A number, which thresholds which points are considered significant. It will also add a red line to the plot to mark significance.
#' @param topLimit  A number, limiting how many points to label. By default, it will label the top 20 points by FDR, if there are more than 20 points.
#' @param max.overlaps The maximum number of overlapping labels, a parameter which is passed to ggrepel when labeling points on the plot. 
#' @param returnDF Boolean flag to determine whether to return the ggplot or a data.frame of the information. If FALSE, it will return a ggplot and not a data.frame. 
#' @param addSignificanceLine Boolean flag to determine whether or not to add a line at the FDR significance threshold. 
#' @param rastr A boolean flag to determine whether to rasterize the pdfs. Requires ggrastr to be installed.
#'
#' @return A ggplot object, showing number of interactions on the y-axis and the rank order of measurements on the x-axis 
#'    (Measurements rank-ordered by the number of interactions they have with the other modality.)
#'
#'  
#'
#'
#' @export
#' @keywords model_results

volcanoPairs <- function(interModel, FDR_threshold = 0.05,  topLimit = 20, max.overlaps = 10, returnDF = FALSE, rastr = FALSE, addSignificanceLine = TRUE){

    FDR <- Estimate <- Pair <- NULL
    interAct <- as.data.frame(SummarizedExperiment::assays(interModel)[['exp2']])
    #Generate range for plotting
    xRange <- range(interAct$Estimate)*0.10
    yRange <- range(-log10(interAct$FDR))*0.10
    #Generate labels for each plot
    interAct$Pair = rownames(interAct)
    #Filter out which pairs to label. First by FDR threshold, then by the topLimit parameter. 
    interAct$Pair[interAct$FDR >= FDR_threshold] = NA
    if(sum(interAct$FDR <= FDR_threshold)){
        topPairs <-  interAct$Pair[order(interAct$FDR)[c(1:topLimit)]]
        interAct$Pair[!interAct$Pair %in% topPairs] = NA
    }

    if(returnDF){
      return(interAct)
    }

    #Generate ggplot.  

    if(rastr){
        if (!requireNamespace("ggrastr", quietly = TRUE)) {
            stop(
                "Package 'ggrastr' is required for volcanoPairs when rastr = TRUE. ",
                "Please install 'ggrastr' or set rastr = FALSE to proceed."
            )
        }

        p1 <- ggplot2::ggplot(interAct, ggplot2::aes(x = Estimate, y = -log10(FDR), label = Pair)) + 
           ggplot2::geom_point_rast() +   ggplot2::coord_cartesian(clip = 'off') + 
                ggplot2::theme(plot.margin = ggplot2::margin(1,1,1.5,1.2, "cm")) +
            ggrepel::geom_text_repel(max.overlaps = max.overlaps, force = 4,
                            max.iter = 100000) + ggplot2::theme_minimal() 

    } else {

        p1 <- ggplot2::ggplot(interAct, ggplot2::aes(x = Estimate, y = -log10(FDR), label = Pair)) + 
           ggplot2::geom_point() +   ggplot2::coord_cartesian(clip = 'off') + 
                ggplot2::theme(plot.margin = ggplot2::margin(1,1,1.5,1.2, "cm")) +
            ggrepel::geom_text_repel(max.overlaps = max.overlaps, force = 4,
                            max.iter = 100000) + ggplot2::theme_minimal() 
    }


  
    if(addSignificanceLine){

      p1 <- p1+ 
            ggplot2::geom_hline(yintercept = -log10(FDR_threshold), color = 'red', alpha = 0.25) 
            
    }
    return(p1)

}



#' @title Extract association model predictions with original data and metadata for plotting
#'
#' @description \code{associatedPredictions} Uses original data and an interaction model to generate group-level interaction predictions for 
#'                  for visualizations. It will generate a data frame of model-adjusted values from SE1, values from SE2, 
#'                  group-level model predictions from values in SE2, and all metadata. The dataframe for each interaction pair
#'                  will be saved within one assay of an SummarizedExperiment and exported.         
#' @param SE1 SummarizedExperiment object of data from measurement-type 1 that was predicted from an interaction model.
#' @param SE2 SummarizedExperiment object for measurement-type 1 
#' @param assay1 - name of the assay to extract and use from SE1
#' @param assay2 - name of the assay to extract and use from SE2
#' @param sampleColumn - A string: the name of the column with Sample names from the metadata. 
#'                  Must be the same from SE1 and SE2. This is used for aligning SE1 and SE2
#' @param interModel The output of an association model (GeneTile, scATAC_Associations, scRNA_Associations, or general_Association)
#' @param firstPair Measurement name from SE1 to look at. The measurement was predicted using the measurement in secondPair, and will be adjusted for all other continuous model factors, except second pair, if adjust = TRUE. 
#' @param secondPair Measurement name from SE2 to look at. This measurement will be the one used to predict the measurment provided to firstPair.  
#' @param interModel A SummarizedExperiment object contained model coefficients and generated from any associative model function (general_Associations, scRNA_Associations, scATAC_Associations, and GeneTile_Associations).
#' @param adjust A boolean, default to TRUE. Determines whether to generate an column with an adjusted measurement. The measurement will be adjusted for all factors that are not the specVariable within the model.
#' @param outputDF A boolea, default is FALSE, as to whether to return a summarized experiment with a seperate assay for each pair, or just a list of data.frames.
#' @return A SummarizedExperiment, where each data.frame is the original data used to associated, with adjusted data (if adjust = TRUE), the model predictions, and all other metadata. 
#'
#' @export
#' @keywords model_results
#'


associatedPredictions <- function(SE1, SE2, assay1, assay2, sampleColumn, interModel, firstPair, secondPair, adjust = TRUE, outputDF =FALSE){

    ## Make sure the SE1 is a ChAIObject and then combine so that sample metadata can be pulled out. 
    se1_type = isChAIObject(SE1, type = 'data', returnType = TRUE)
    se2_type = isChAIObject(SE2, type = 'data', returnType = TRUE)
    model_type = isChAIObject(interModel, type = 'model', returnType = TRUE)
    
    logLink <- model_type == 'scRNA_Assocations'
    
    #Check whether sampleColumn is present in both SE objects. 
    if(!(sampleColumn %in% colnames(SE1@colData) & sampleColumn %in% colnames(SE2@colData))){
      stop('sampleColumn not found in SE1 and/or SE2.')
    }
    
    ##Check data type to make sure SE1 and SE2 are in the right order. 
    if(model_type == 'GeneTile' & c(se2_type != 'scRNA' | se1_type != 'scATAC')){
      stop('A GeneTile model provided. SE1 should be a MOCHA SampleTile Object containing pseudobulking scATAC data, and 
                  SE2 should be ChAI scRNA object with normalized pseudobulk RNA. SE1 and/or SE2 are not matching expectations.')
      
    }else if(model_type == 'scRNA_Associations' & se1_type != 'scATAC'){
      stop('A scRNA Association model provided. SE1 should be ChAI scRNA object with normalized pseudobulk RNA and it is not.')
      
    }else if(model_type == 'scATAC_Associations' & se1_type != 'scATAC'){
      stop('A scATAC Association model provided. SE1 should be a MOCHA SampleTile Object containing pseudobulking scATAC data, but it is not.')
      
    }else if(model_type == 'ChromVAR_Associations' & se1_type != 'ChromVAR'){
      stop('A ChromVAR Association model provided. SE1 should be a ChAI ChromVAR object, generated by makeChromVAR,  but it is not.')
      
    }else if(model_type == 'General_Associations' & c(se2_type != 'General' | se1_type != 'General')){
      
      stop('A General Association model provided. Both SE1 and SE2 should be general ChAI data objects, from importGeneralModality, but they are not.')
      
    }else if(!grepl('GeneTile|Association', model_type)){
      
      stop("interModel does not appear to be an associative model. Please check that the variable interModel is correct")
      
    }
    
    if(model_type == 'General_Associations'){
      warning('associatedPredictions does notcurrently support General Association models with non-Gaussian error distributions.')
    }
    
    #Check whether samples align.     
    if(all(!SE1@colData[,sampleColumn] %in% SE2@colData[,sampleColumn]) |
            all(!SE2@colData[,sampleColumn] %in% SE1@colData[,sampleColumn])){

        stop(stringr::str_interp('samples names in ${sampleColumn} are not the same. Please ensure that sample names match in SE1 and SE2 via ${sampleColumn}.'))

    } else if(!all(SE2@colData[,sampleColumn] %in% SE1@colData[,sampleColumn]) |
        !all(SE1@colData[,sampleColumn] %in% SE2@colData[,sampleColumn])) {
        
        partialMatch <- intersect(SE2@colData[,sampleColumn] , 
                                    SE1@colData[,sampleColumn])

        warning(stringr::str_interp('Some samples names in ${sampleColumn} are not the same between modalities. Non-matching names will be dropped'))
        dropped = SE2@colData[,sampleColumn] %in% SE1@colData[,sampleColumn]
        generalDropped = SE1@colData[,sampleColumn] %in% SE2@colData[,sampleColumn]
        warning(stringr::str_interp(' ${sum(!dropped)} and ${sum(!generalDropped)} sampled dropped from SE2 and SE1, respectively'))

        if(se1_type == 'scATAC'){

            SE1 <- subsetMOCHAObject(SE1, subsetBy = sampleColumn, 
                                    groupList = partialMatch)

        }else{

            SE1 <-  subsetChAI(SE1, subsetBy = sampleColumn, 
                                groupList = partialMatch)

        }

        SE2 <-  subsetChAI(SE2, subsetBy = sampleColumn, 
                                groupList = partialMatch)
        
    }


    if(!all(SE1@colData[,sampleColumn] == SE2@colData[,sampleColumn])){

        warning('Reording rnaSE sample names to match atacSE.')
        SE2 <- SE2[,match(SE1@colData[,sampleColumn], SE2@colData[,sampleColumn])]

    }

    #Test whether the pair is found within their assays
    if(!all(firstPair %in% rownames(SE1))){
      stop('firstPair not found within SE1. Please read documentation.')
    }
    if(all(secondPair %in% rownames(SE1))){
      stop('Features from secondPair were found in SE1. Please ')
    }else if(!all(secondPair %in% rownames(SE2))){
      stop('secondPair not found within SE2. Please read documentation.')
    }
    
    if(length(secondPair) != length(firstPair)){
      stop('firstPair and secondPair should be the same length. This function does not find all combinations of firstPair and secondPair, but only matched pairs (i.e. A1, B2, B3, not A1, A2, A3')
    }
    
    pair = as.list(paste(firstPair, secondPair, sep ='_'))
    #Test whether the pair was actually modeled
    if(!all(pair %in% rownames(interModel))){
      stop('Pair of interacting features not found within interModel. This pair of features was not actually modeled.')
    }
    
    if(any(c(colnames(SummarizedExperiment::colData(SE1)), 
             colnames(SummarizedExperiment::colData(SE2))) %in% c('exp2'))){
      
      stop('metadata of SE1 and/or SE2 contains a column that contains the name exp2.',
           'exp2 are hardcoded to represent the data from SE1 and SE2, not the metadata.
          Please remove these and try again.')
      
    }
    
    if(!assay1 %in% names(SummarizedExperiment::assays(SE1))){
      stop('assay1 not found in SE1')
    }
    
    if(!assay2 %in% names(SummarizedExperiment::assays(SE2))){
      stop('assay2 not found in SE2')
    }
    
    
    ## Pull out data and metadata.
    if(se1_type %in% c('scRNA', 'ChromVAR')){
      
      mat1 <- SummarizedExperiment::assays(SE1[firstPair,])[[assay1]] 
      mat1[is.na(mat1)] = 0
      SE12 = flattenChAI(SE1, cellPopulations = assay1, metadataT= TRUE)
      metaData <-  SummarizedExperiment::colData(SE12)
      
    }else if(se1_type %in% c('scATAC')){
      
      mat1 <- SummarizedExperiment::assays(SE1[firstPair,])[[assay1]] 
      mat1[is.na(mat1)] = 0
      subSTM <- MOCHA::subsetMOCHAObject(SE1, subsetBy = 'celltype', groupList = assay1)
      SE12 = MOCHA::combineSampleTileMatrix(subSTM)
      metaData <-  SummarizedExperiment::colData(SE12)
      
    }else{
      metaData <-  SummarizedExperiment::colData(SE1)
      mat1 <- SummarizedExperiment::assays(SE1[firstPair,])[[assay1]] 
      
      mat1[is.na(mat1)] = 0
    }
    mat2 <- SummarizedExperiment::assays(SE2[secondPair,])[[assay2]] 

    #Check data type to determine transformations. 
    if(model_type %in% c('scATAC_Associations','GeneTile')){
        mat1 <- log2(mat1 + 1)
    }
    if(model_type == 'GeneTile'){
        mat2 <- log2(mat2 + 1)
     }

    assayList <- SummarizedExperiment::assays(interModel)
    allVariables <- names(assayList)[!(names(assayList) %in% c('exp2') | grepl('Intercept', names(assayList)))]
    if(any(grepl('ZI_', allVariables))){
    
        warning('Zero-inflated variables found in modelSE. This function works with continuous fits. Zero inflated model estimates will be dropped')
        allVariables = allVariables[!grepl('ZI_', allVariables)]
        
    }
    numericVariables <- names(assayList)[names(assayList) %in% colnames(metaData)]

    remainingVariables <- allVariables[!allVariables %in% numericVariables & !grepl('ZI', allVariables)]
       
    ## Log transform the FragNumbers & CellCounts so as to stabilize the model. But only if FragNumber is in the model. Same for CellCounts.
    if(any(colnames(metaData ) %in% c('FragNumber'))){
        metaData$rawFragNumber = metaData$FragNumber
        metaData$FragNumber <- log10(metaData $FragNumber)
    }
    if(any(colnames(metaData ) %in% c('CellCounts'))){
        metaData$rawCellCounts = metaData$CellCounts
        metaData$CellCounts <- log10(metaData $CellCounts)
    }

    #Create a metadata column for each categorical variable, one-hot encoding them. 
    if(length(remainingVariables) > 1){

        nextVariables <- lapply(remainingVariables, function(x) matchCategorical(metaData, x))
        newMetaData = do.call('cbind', nextVariables)
        colnames(newMetaData) <- remainingVariables
        metaData = cbind(metaData, newMetaData)

    }else if(length(remainingVariables) > 0){
        newMetaData <- data.frame(newVar = matchCategorical(metaData, remainingVariables))
        colnames(newMetaData) <- remainingVariables
        metaData = cbind(metaData, newMetaData)
    }


    subMeta <- metaData[, colnames(metaData) %in% allVariables, drop = FALSE]
                         
    allPredictions <- pbapply::pblapply(cl = NULL, seq_along(pair), function(x){
        
        modelVals <- getModelValues(interModel, pair[[x]])[,'Estimate', drop=FALSE]
        
        bothData <- data.frame('exp1' =  unlist(mat1[firstPair[[x]],]), 'exp2' = unlist(mat2[secondPair[[x]],]))
        rownames(bothData) = colnames(mat1)
        
        
        if(adjust & length(allVariables) > 0){

            bothData$orig_exp1 = bothData$exp1
            allAdjusts = as.matrix(subMeta[,allVariables] ) %*% as.matrix(modelVals[allVariables,,drop=FALSE])
            predictVariables <- as.matrix(bothData[, 'exp2', drop = FALSE]) %*% 
              as.matrix(modelVals['exp2', , drop = FALSE])
            
            if(logLink){
                logExp1 <- log(bothData$exp1) - rowSums(allAdjusts)
                bothData$exp1 = exp(logExp1)
            } else{
                bothData$exp1 = bothData$exp1 - rowSums(allAdjusts)  
            }
            
             if(any(grepl('Intercept', rownames(modelVals))) & !logLink){
                 
                bothData$Prediction = modelVals[grepl('Intercept', rownames(modelVals)),] + 
                                        unlist(as.list(predictVariables))
                 
            }else if(any(grepl('Intercept', rownames(modelVals)))){
                
                bothData$Prediction = exp(modelVals[grepl('Intercept', rownames(modelVals)),] + 
                                        unlist(as.list(predictVariables)))
                
            } else{
                 
                bothData$Prediction = exp(unlist(as.list(predictVariables)))
                 
            }
                 
            bothData$Prediction = modelVals[grepl('Intercept',rownames(modelVals)),] + modelVals['exp2',]*bothData$exp2
        
        }else if(length(allVariables) > 0){
    
            numericCalls <- sum(unlist(lapply(numericVariables, function(x) {
                                modelVals[x,]*mean(subMeta[,x])
                })))
            
            bothData$Prediction = modelVals['Intercept',] + modelVals['exp2',]*bothData$exp2 + numericCalls
            
        }else{
            
            bothData$Prediction = modelVals['Intercept',] + modelVals['exp2',]*bothData$exp2
            
        }
            
        if(logLink){
            
                newData$Prediction = exp(newData$Prediction)
                
        }
            
        bothData <- cbind(bothData, metaData)
        bothData
        
        #Check data type to determine transformations. 
        if(model_type %in% c('scATAC_Associations','GeneTile') & adjust){
            bothData$exp1 = 2^(bothData$exp1)-1
            bothData$orig_exp1 = 2^(bothData$orig_exp1)-1
        }else if(interModel@metadata$Type %in% c('scATAC_Associations','GeneTile')){
            bothData$exp1 = 2^(bothData$exp1)-1
        }
            
        if(model_type %in% 'GeneTile'){
            bothData$exp2 = 2^(bothData$exp1)-1

        }
        bothData
        
    })
        
    names(allPredictions) = pair
                                
    if(length(allPredictions) == 1 & outputDF){
    
        return(allPredictions[[1]])
        
    }else if(outputDF){
    
        return(allPredictions)
    }

    exp = SummarizedExperiment::SummarizedExperiment(allPredictions, metadata = interModel@metadata)
   
    return(exp)
}

#' @title Export Gene-tile associations for visualization or browsing in IGV-viewer
#'
#' @description \code{exportGeneTileLinks} Takes the output of GeneTile_Associations and exports a file in the bedpe format for easy browsing of gene-tile links.
#' @param GeneTileObj SummarizedExperiment objectoutput of GeneTile_Associations.
#' @param fileName String describing the name of the file to write to disk. 
#' @param returnDF Boolean. Default is FALSE. Determines whether to return all gene-tile links after writing out a file. 
#' @param FDR_threshold A number between 0 and 1. Determines threshold for calling a gene-tile association as significant. Passed to getEstimates
#' @param filterFactors A vector of characters matching variables within the GeneTile model. This is passed to getEstimates, which will remove all GeneTile links that are also significantly associated with any of these factors. 
#' @param backgroundThreshold A number between 0 and 1. Determines threshold for calling a background variable as significant. Passed to getEstimates for filtering out regions. 
#'
#' @export
#' @keywords networks
#'

exportGeneTileLinks <- function(GeneTileObj, fileName,
                                returnDF = FALSE,
                                FDR_threshold = 0.1,
                                filterFactors = NULL, 
                                backgroundThreshold =0.1){
    
        modelType = isChAIObject(GeneTileObj, type = 'model', returnType = TRUE)
    
        if(modelType !=  'GeneTile'){
                
                stop('GeneTileObj is not actually a GeneTile Association from ChAI.')
                
        }
    
        if(!methods::is(fileName, 'character')){
        
            stop('fileName is not a character')
            
        }

    
        geneTile1 <- getEstimates(GeneTileObj, factor = 'exp2', FDR_threshold= FDR_threshold,
                                 filterFactors = filterFactors, backgroundThreshold = backgroundThreshold)
        geneTile1$Type = 'Continuous'
        geneTile2 <- getEstimates(GeneTileObj, factor = 'ZI_exp2', FDR_threshold= FDR_threshold,
                                 filterFactors = filterFactors, backgroundThreshold = backgroundThreshold)
        geneTile2$Type = 'ZeroInflated'
    
        allGeneTiles = rbind(geneTile1, geneTile2)
    
        allRowData <- SummarizedExperiment::rowData(GeneTileObj)
    
        #Filter row data to just gene tiles of interest.
        allRowData <- allRowData[allGeneTiles$Obj,]
    
        allLinks = cbind(
                    data.frame(chr1 = gsub(":.*", "", allRowData$Promoters),
                            start1 = gsub("-.*","", gsub(".*:", "", allRowData$Promoters)),
                            end1 = gsub(".*-", "", allRowData$Promoters),
                            chr2 = gsub(":.*", "", allRowData$Obj1),
                            start2 = gsub("-.*","", gsub(".*:", "", allRowData$Obj1)),
                            end2 = gsub(".*-", "", allRowData$Obj1),
                            name = rownames(allRowData),
                            score = -log10(allGeneTiles$FDR),
                            strand1 = rep("*", length(allRowData$Promoters)),
                            strand2 = rep("*", length(allRowData$Promoters))
                         ),
                      allGeneTiles
                       )
    
        if(any('tileType' == rownames(allRowData))){
        
            allLinks$tileType = allRowData$tileType
            allLinks$tileGene = allRowData$Gene
        }
    
       utils::write.table(allLinks, fileName, 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep ='\t')
    
        if(returnDF){
        
            return(allLinks)
           
        }
 
}


#' @title Get a data.frame with of model coefficients 
#'
#' @description \code{getSpecificPair} Takes the output of associatedPredictions and returns a data.frame for a single combination of objects.        
#' @param associationPredictions SummarizedExperiment object output from associatedPredictions. 
#' @param firstPair A single measurement name from SE1
#' @param secondPair A single measurement name from SE2 
#' @return A data.frame of original data from SE1, adjusted data (if adjust = TRUE), model predictions, and all other metadata. 
#'
#' @export
#' @keywords model_results


getSpecificPair <- function(associationPredictions, firstPair, secondPair){

    association_type = isChAIObject(associationPredictions, type = 'model', returnType = TRUE)

    allFirst <- gsub("_.*", "", names(SummarizedExperiment::assays(associationPredictions)))
    allSecond <- gsub(".*_", "", names(SummarizedExperiment::assays(associationPredictions)))

    if(any(length(firstPair) + length(secondPair) != 2)){
        stop('firstPair and secondPair should only contain one string. This function is for getting a data.frame of one specific predicted association.')
    }

    if(any(!firstPair %in% allFirst)){
        stop('The firstPair value is not found in the associatedPredictions object. Are you sure it was tested?')
    }

    if(any(!secondPair %in% allSecond)){
        stop('The secondPair value is not found in the associatedPredictions object. Are you sure it was tested?')
    }

    

    pair = paste(firstPair, secondPair, sep = "_")

    if(!pair %in% names(SummarizedExperiment::assays(associationPredictions))){
        stop(stringr::str_interp('The combination of ${firstPair} and ${secondPair} are not in the associatedPredictions. Please regenerate the associatedPredictions object with these two measurements, or try a different combination.'))
    }

    mat1 <- as.data.frame(SummarizedExperiment::assays(associationPredictions)[[pair]])

    if(associationPredictions@metadata$Type %in% c('GeneTile')){

        colnames(mat1)[1] = 'Tiles'
        colnames(mat1)[2] = 'Genes'
        mat1$TileName = firstPair
        mat1$GeneName = secondPair

    } else if(associationPredictions@metadata$Type %in% c('scATAC_Associations')){

        colnames(mat1)[1] = 'Tiles'
        colnames(mat1)[2] = 'Modality2'
        mat1$TileName = firstPair
        mat1$MeasurementName = secondPair

    } else if(associationPredictions@metadata$Type %in% c('scRNA_Associations')){

        colnames(mat1)[1] = 'Genes'
        colnames(mat1)[2] = 'Modality2'
        mat1$GeneName = firstPair
        mat1$MeasurementName = secondPair

    }  else if(associationPredictions@metadata$Type %in% c('ChromVAR_Associations')){

        colnames(mat1)[1] = 'TFs'
        colnames(mat1)[2] = 'Modality2'
        mat1$GeneName = firstPair
        mat1$MeasurementName = secondPair

    } else if(associationPredictions@metadata$Type %in% c('General_Associations')){

        colnames(mat1)[1] = 'Modality1'
        colnames(mat1)[2] = 'Modality2'
        mat1$FirstMeasurementName = firstPair
        mat1$SecondMeasurementName = secondPair
        
    }else{

        stop('association type not recognized. Please check whether this was generated by a ChAI association function.')

    }

    return(mat1)

}


#' @title matchCategorical 
#'
#' @description \code{matchCategorical} Identifies categorical variables from the metadata and one-hot encodes them. 
#' @param  metaDataDF
#' @return A one-hot encoding for the categorical variable given, based on the metadata. There will be one column for all options. 
#'
#'
#' @noRd
matchCategorical <- function(metaDataDF, variable, spacer = '.'){
    
    metadata2 = metaDataDF
    colnames(metadata2)= paste("ZI_", colnames(metadata2), sep = "_")
    metadata3 = cbind(metadata2, metaDataDF)
    whichMatch = which(unlist(lapply(colnames(metaDataDF), function(x) grepl(x, variable))))
    if(length(whichMatch) > 1){
       
       #Iterate over all possible matching columns, and see if any combinations of column name and 
       #values match the variable. 
       allMatches <- lapply(whichMatch, function(x){

            values = gsub(" ",spacer,unlist(metaDataDF[, whichMatch]))
            valuesList <- paste(colnames(metaDataDF)[whichMatch], values, sep ='')
            if(!any(valuesList == variable)){

               NA

            }else{

              valuesList == variable
            } 
           
        })
        #See if multiple columns match. If multiple match, then see if they are all the same or different.
        # if different, throw an error. If the same or only one match, then use it. 
        whichAllMatch <- unlist(lapply(allMatches, function(x) all(is.na(x))))
        if(sum(!whichAllMatch) > 1){

           combMatch <- do.call('rbind', allMatches[!whichAllMatch])
           if(any(unlist(lengths(apply(combMatch, 2, unique))) > 1)){
           
                stop(paste('The categorical variable ', variable,
                    ' could match with multiple different metadata columns. Please clean up metadata.', sep=''))
        
           }else{
               
                return(allMatches[[1]])
                
            }

        }else if(sum(!whichAllMatch) == 1){
            
          return(allMatches[[!whichAllMatch]])
            
        }else{
            
          stop(paste('The categorical variable ', variable, ' could not be matched to metadata', sep=''))
            
        }


   }else if(length(whichMatch) == 1){
       
       #If only one column matches, then pull that one out. 
        values = gsub(" ",spacer,unlist(metaDataDF[, whichMatch]))
        valuesList <- paste(colnames(metaDataDF)[whichMatch], values, sep ='')
       if(!any(valuesList == variable)){
            
           stop(paste('The categorical variable ', variable, ' could not be matched to metadata', sep=''))
           
       }else{
           
           outputList <-  valuesList == variable
           return(outputList)
       }
       
       
   }else{
       
       stop('Variable not found in metadata. Please ensure all variables for the model are within the metadata of the summarized experiment.')
       
   }

}
