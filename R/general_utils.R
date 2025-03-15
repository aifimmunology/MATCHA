#' @title Create a ChAI General Modality object from any feature count matrix and sample metadata.
#'
#' @description \code{importGeneralModality} Accepts normalized data from any general modality, associated metadata, and generates a ChAI-compatible SummarizedExperiment object. 
#' @param DataMatrix A data.frame, data.table, or matrix, where rows are measurements and columns are samples. Column names should match the names found within the 'Sample' column of the MetaData you provide. 
#'  Alternatively, DataMatrix can also be named list of matrices, where all matrices are the same size, with the same row names and column names. These could be, for example, a list of gene expression across samples for multiple sorted cell types.
#' @param MetaData Associated metadata for each sample within the DataMatrix. Must contain a column named 'Sample' that has all the columns names within the DataMatrix
#' @param sampleColumn String. Name of the column that has Sample-level metadata. Default is Sample. Will overwrite any column that has Sample already present. 
#' @return A SummarizedExperiment with all your data inerwoven with the MetaData
#'
#'
#' @keywords data_import
#'
#' @export


importGeneralModality <- function(DataMatrix, MetaData, sampleColumn = 'Sample'){
    
    
    if(!methods::is(MetaData, 'data.frame')){
        stop('MetaData is not a data.frame. MetaDatamust be a data.frame object.')
    }

    if( methods::is(DataMatrix, 'data.frame') |  methods::is(DataMatrix, 'data.table') |  methods::is(DataMatrix, 'matrix')){
        newMatrix <- as.data.frame(DataMatrix)
        rownames(newMatrix) <- rownames(DataMatrix)
        specSamples <- colnames(DataMatrix)

        if(!'Sample' %in% colnames(MetaData) & sampleColumn %in% colnames(MetaData)){
            MetaData$Sample = unlist(MetaData[,sampleColumn])
        }else if(!sampleColumn %in% colnames(MetaData)){

            stop('sampleColumn does not exist in MetaData')

        }else if('Sample' %in% colnames(MetaData)){
        
            if(!all(MetaData[,sampleColumn] == MetaData[,'Sample'])){
            warning(stringr::str_interp("${sampleColumn} and Sample column in the MetaData data do not match. Over-writing MetaData column named 'Sample' since sampleColumn = 'Sample'"))
                MetaData$Sample = MetaData[,sampleColumn]
            }
            
        }else{
        
            MetaData$Sample = MetaData[,sampleColumn]
            
        }
            
        #Reorder metadata and matrix
        newMeta <- MetaData[order(MetaData$Sample),]
        newMatrix <- newMatrix[,order(specSamples)]
        rownames(newMeta) = newMeta$Sample 
        specSamples <- specSamples[order(specSamples)]
                                         
        if(!all(MetaData$Sample %in% specSamples) | !all(specSamples %in% MetaData$Sample)){
            stop('DataMatrix and MetaData do not align. Please ensure all column names for dataMatrix can be found within the Sample column of MetaData.')
        }else{
            
            newSE <- SummarizedExperiment::SummarizedExperiment(list('General' = newMatrix),
                            colData = newMeta, 
                            metadata = list('Type' = 'General', 
                                'History' = paste("importGeneralModality", 
                                                  utils::packageVersion("ChAI"))
                                           )
                        )

        }

    }else if(methods::is(DataMatrix, 'list')){


        specSamples <- unique(lapply(DataMatrix, colnames))
        allDims <- do.call('rbind',lapply(DataMatrix,dim))
        #Are all the row lenghts and column lengths the same?
        if(length(unique(allDims[,1])) != 1 | length(unique(allDims[,2])) != 1  ){

            stop('You have provided a list of matrices, but they are not all the same dimensions. Please make sure they all the same dimensions.')

        }
            
        #Set the SampleColumn in MetaData so that the data and metadata are aligned.
        if(!'Sample' %in% colnames(MetaData) & sampleColumn %in% colnames(MetaData)){
            MetaData$Sample = unlist(MetaData[,sampleColumn])
        }else if(!sampleColumn %in% colnames(MetaData)){

            stop('sampleColumn does not exist in MetaData')

        }else if('Sample' %in% colnames(MetaData)){
        
            if(!all(MetaData[,sampleColumn] == MetaData[,'Sample'])){
            warning(stringr::str_interp("${sampleColumn} and Sample column in the MetaData data do not match. Over-writing MetaData column named 'Sample' since sampleColumn = 'Sample'"))
                MetaData$Sample = MetaData[,sampleColumn]
            }
            
        }else{
        
            MetaData$Sample = MetaData[,sampleColumn]
            
        }
            
        newMeta <- MetaData[order(MetaData$Sample),]
        rownames(newMeta) = newMeta$Sample
        newMatrix <- lapply(DataMatrix, function(x) x[,order(specSamples)])
        specSamples <- specSamples[order(specSamples)]
                                         
        #Check to make sure the same sample names are present and then make the Summarized Experiment.
        if(!all(MetaData$Sample %in% specSamples) | !all(specSamples %in% MetaData$Sample)){
           stop('DataMatrix and MetaData do not align. Please ensure all column names for dataMatrix can be found within the Sample column of MetaData.')
        }else{
            
            newSE <- SummarizedExperiment::SummarizedExperiment(newMatrix,
                            colData = newMeta, 
                            metadata = list('Type' = 'General', 
                                     'History' = list(paste("importGeneralModality", utils::packageVersion("ChAI")))
                                     )
                            )

        }
    }else{
        stop('DataMatrix not recognized. Must be either a data.frame, data.table, or matrix. or a list of those.')
    }

    return(newSE)

}



#' @title get significant estimates for a given variable from any Model or Association Model
#'
#' @description \code{getEstimates} Pull out a data.frame of filtered estimates with associated STD deviation, p-values, and 
#'              FDR for a single factor from the output SummarizedExperiment object of any modeling function. 
#' @param modelSE A SummarizedExperiment of model associations
#' @param factor A string. Default is exp2, which will pull out associations. If you used a factor in both the Zero-inflated and non-zero inflated models, you must pull each one out individualy. 
#' @param FDR_threshold Threshold of determining significance.Default is 0.1 
#' @param filterFactors A optional list of strings, describing the factors you do not to select for. If a given measurement is below the backgroundThreshold for that factor, it will be removed.
#' @param backgroundThreshold An optional background threshold that determines how stringently you want to remove measurements if they related to any of the filterFactors, based on their FDR. Standard threshold is 0.1

#'
#' @return data frame of coefficient values for a given factor, filtered by an FDR threshold and by the lack of relationship to other factors, if provided. 
#'
#'
#'
#' @examples
#' \dontrun{
#'   age_df <- getEstimates(model_scATAC_output, factor = 'Chr1:500-999', 
#'                                FDR_threshold = 0.1, backgroundThresold = 0.2)
#' }
#'
#' @export
#' @keywords model_results

getEstimates <- function(modelSE, factor = 'exp2', FDR_threshold = 0.1, 
                            filterFactors = NULL, backgroundThreshold = 0.1){

    #Make sure FDR_threshold and backgroundThreshold are numeric between 0 and 1.
    if(!is.numeric(FDR_threshold)){
        stop('FDR_threshold must be numeric.')
    }else if(FDR_threshold > 1 | FDR_threshold < 0){
         stop('FDR_threshold must be numeric between 0 and 1.')
    }
        
    if(!is.numeric(backgroundThreshold)){
        stop('backgroundThreshold must be numeric.')
    }else if(backgroundThreshold > 1 | backgroundThreshold < 0){
         stop('backgroundThreshold must be numeric between 0 and 1.')
    }

    model_type <- isChAIObject(modelSE, type= 'model', returnType = TRUE)

    if(!any(names(SummarizedExperiment::assays(modelSE)) %in% factor)){
        stop(stringr::str_interp('${factor} was not modeled within this modelSE provided.'))
    }

    mat1 <- SummarizedExperiment::assays(modelSE)[[factor]]
    mat1$Obj <- rownames(mat1)

    if(model_type == 'scATAC_Model'){

        colnames(mat1)[colnames(mat1) == 'Obj'] = 'Tiles'

    }else if(model_type %in% c('scRNA_Model')){

        colnames(mat1)[colnames(mat1) == 'Obj'] = 'Genes'

    } else if(model_type %in% c('ChromVAR_Model')){

        colnames(mat1)[colnames(mat1) == 'Obj'] = 'TFs'

    } else if(model_type %in% c('General_Model')){

        colnames(mat1)[colnames(mat1) == 'Obj'] = 'Measurement'

    } else if(model_type %in% c('GeneTile')){

        mat1$Tiles = gsub("_.*","", mat1$Obj)
        mat1$Genes = gsub(".*_","", mat1$Obj)

    } else if(model_type %in% c('scATAC_Associations')){

        mat1$Tiles = gsub("_.*","", mat1$Obj)
        mat1$General = gsub(".*_","", mat1$Obj)

    } else if(model_type %in% c('scRNA_Associations')){

        mat1$Genes = gsub("_.*","", mat1$Obj)
        mat1$General = gsub(".*_","", mat1$Obj)

    }else if(model_type %in% c('ChromVAR_Associations')){

        mat1$TFs = gsub("_.*","", mat1$Obj)
        mat1$General = gsub(".*_","", mat1$Obj)

    } else if(model_type %in% c('General_Associations')){

        mat1$Modality1 = gsub("_.*","", mat1$Obj)
        mat1$Modality2 = gsub(".*_","", mat1$Obj)
        
    }else{

        stop('association type not recognized. Please check whether this was generated by a ChAI association function.')

    }

    if(!is.null(filterFactors)){

        if(!is.numeric(backgroundThreshold)){
            stop('backgroundThreshold must be numeric.')
        }

        if(!all(filterFactors %in% names(SummarizedExperiment::assays(modelSE)))){
            notFound <- filterFactors[which(filterFactors %in% names(SummarizedExperiment::assays(modelSE)))]

            stop(stringr::str_interp('The filterFactor ${notFound} was not modeled within this modelSE provided.'))
        }
    
        measurementsToRemove <- rowSums(do.call('cbind', lapply(filterFactors, function(XX){

            backMat <- as.data.frame(SummarizedExperiment::assays(modelSE)[[XX]])
            backMat$FDR < backgroundThreshold & !is.na(backMat$FDR)

        })))
        
        #Now filter by FDR
        mat1 <- mat1[mat1$FDR < FDR_threshold & !is.na(mat1$FDR) & measurementsToRemove == 0,]

    }else{
       
        #Now filter by FDR
        mat1 <- mat1[mat1$FDR < FDR_threshold & !is.na(mat1$FDR),]
        
    }


    
    return(mat1)
}


#' @title Get a list of factors in a model
#'
#' @description \code{getModelFactors} Takes a model object and returns all the names of the factors used. 
#' @param modelObj A model object output from any ChAI modeling function (whether on a single modality or the association between two modalities)
#' @return A list of strings. 
#'
#'
#'
#' @examples
#' \dontrun{
#'   getModelFactors(modelObj)
#' }
#'
#' @export
#' @keywords model_results
#'

getModelFactors <- function(modelObj){

    return(names(SummarizedExperiment::assays(modelObj)))

}


#' @title Get a list of all support types of ChAI models
#'
#' @description \code{getModelTypes} Returns a list of all potential model types, or checks what type a given model is. 
#' @param optional string
#'
#' @return A string
#'
#'
#'
#' @examples
#' \dontrun{
#'   age_df <- getModelTypes()
#' }
#'
#' @export
#' @keywords model_results
#'

getModelTypes <- function(optional = NULL){

    if(is.null(optional)){

        modelList <- c('General_Associations','scATAC_Associations', 'scRNA_Associations','ChromVAR_Associations', 'GeneTile','General_Model','ChromVAR_Model', 'scATAC_Model','scRNA_Model')
        return(modelList)

    }
    if(isChAIObject(optional, type = 'model')){

        return(optional@metadata$Type)
    }

}