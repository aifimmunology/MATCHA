#' @title Generate a ChAI ChromVAR Object from MOCHA's Tile-by-Sample Matrix Object
#'
#' @description \code{makeChromVAR} Runs ChromVAR on MOCHA's Tile-by-Sample Matrix object. Slow funciton. 
#'
#' @param atacSE A MOCHA Tile-by-Sample Object (SummarizedExperiment)generated from getSampleTileMatrix within \code{\link[MOCHA]{MOCHA}}. 
#' @param cellPopulation Names of cell types to analyze. Must match assay names in the SummarizedExperiment object. 
#'  Alternative, if you want to run ChromVAR across all celltypes, you can provide the output of combineSampleTileMatrix, and set this parameter to 'counts'.
#' @param motifName Name of metadata slot that has motif information for analysis.
#' @param withinCellType Boolean. Default is FALSE, in which case accessibility across celltypes will be merged into one matrix for ChromVAR, rather than running ChromVAR seperately on the matrix for each cell type when withinCellType = TRUE.
#' @param exportRaw Boolean. Default is FALSE, and will export a ChAI object with ChromVAR deviations reformated for modeling. If set to true, it will either export the original chromVARDeviations object for analysis, either as a list by cell type if withinCellType is TRUE or as a single object if withinCellType is TRUE. This chromVARDeviations object can be later reformatted for ChAI using reformatChromVAR.
#' @param numCores Default is 1. Uses chromVAR's standard parallelization, which has memory leak errors at times. Use at your own risk. 
#' @param verbose Boolean.
#' @return A named list of ChromVAR objects. 
#' 
#' @keywords data_import
#'
#' @export

makeChromVAR <- function(atacSE, motifName,
                      cellPopulation = NULL,
                      withinCellType = FALSE,
                      exportRaw = FALSE,
                    numCores = 1,
                      verbose = TRUE) {

    if (!requireNamespace("chromVAR", quietly = TRUE)) {
        stop(
        "Package 'chromVAR' is required for makeChromVAR. ",
        "Please install 'chromVAR' to proceed."
        )
    }
    
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
        stop(
        "Package 'BiocParallel' is required for makeChromVAR. ",
        "Please install 'BiocParallel' to proceed."
        )
    }

    if(is.null(cellPopulation)){

        cellPopulation = names(SummarizedExperiment::assays(atacSE))

    }

    if(is.null(cellPopulation)){

        cellPopulation = names(SummarizedExperiment::assays(atacSE))

    }

    if (
        any(!cellPopulation %in% names(SummarizedExperiment::assays(atacSE)))
    ) {

        stop("cellPopulation was not found within atacSE.")

    } else if (withinCellType) {
   
        #Generate a list of combined objects
        newObj <- lapply(cellPopulation, function(x){
             MOCHA::combineSampleTileMatrix(MOCHA::subsetMOCHAObject(atacSE, subsetBy = 'celltype', groupList = x, subsetPeaks = TRUE))
        })
        names(newObj) <- cellPopulation

    } else if(all(cellPopulation == 'counts')){
        newObj <- atacSE
    }else{
        newObj <- MOCHA::combineSampleTileMatrix(MOCHA::subsetMOCHAObject(atacSE, 
                                                                          subsetBy = 'celltype', groupList = cellPopulation, subsetPeaks = TRUE))
    }
    
    genome <- S4Vectors::metadata(atacSE)$Genome
    genome <- BSgenome::getBSgenome(genome)

    if(numCores > 1){
        
        BiocParallel::register(BiocParallel::SerialParam())
        
    }else if(Sys.info()[['sysname']] != 'Windows'){
        
        BiocParallel::register(BiocParallel::MulticoreParam(numCores, progressbar = TRUE))
        
    }else{
        
        BiocParallel::register(BiocParallel::SnowParam(workers = numCores, type = "SOCK"))
        
    }

    
    #Either iterate ovewr a list of newObj, or just directly on newObj to generate ChromVAR
    if(any(tolower(class(newObj)) %in% 'list')){

        chromVAROut <- lapply(cellPopulation, function(XX){
            if(verbose){ message('Analyzing ', XX)}
            newObj[[XX]] <- addGCBias_ChAI(newObj[[XX]], genome, verbose)
            anno_ix <- chromVAR::getAnnotations(newObj[[XX]]@metadata[[motifName]], 
                        rowRanges = SummarizedExperiment::rowRanges(newObj[[XX]]))
            chromVAR::computeDeviations(object =newObj[[XX]], 
                annotations = anno_ix)
        })
         names(chromVAROut) <- cellPopulation
        
        if(exportRaw){
            return(chromVAROut)
        }

    }else{

        if(verbose){ message('Analyzing ', cellPopulation)}
        newObj <- addGCBias_ChAI(newObj, genome, verbose)

        anno_ix <- chromVAR::getAnnotations(newObj@metadata[[motifName]], 
                        rowRanges = SummarizedExperiment::rowRanges(newObj))
        chromVAROut <- chromVAR::computeDeviations(object = newObj, 
                    annotations = anno_ix)
        
        if(exportRaw){
            return(chromVAROut)
        }
    
    }

    newOut_Dev <- reformatChromVAR(chromVAROut, selectDev =TRUE)
    newOut_Z <- reformatChromVAR(chromVAROut, selectDev =FALSE)

    newOut <- list('Z_Score' = newOut_Z, 'Deviations' = newOut_Dev)

    BiocParallel::register(BiocParallel::SerialParam())


    return(newOut)

}




#' @title \code{reformatChromVARList}
#'
#' @description \code{reformatChromVARList} pseodubulks a Seurat object by sample and cell type into a SummarizedExperiment, similar to MOCHA
#'
#' @param atacSE A SummarizedExperiment object generated from
#'   getSampleTileMatrix or combineSampleTileMatrix
#' @param genome BS Genome file to use as reference for GC bias scoring. 
#' @param verbose Boolean. verbose flag. Default is FALSE. 
#' @return A chromVAR object with GC bias
#'
#' @noRd

addGCBias_ChAI <- function(obj1, genome, verbose = FALSE){

    if (!requireNamespace("chromVAR", quietly = TRUE)) {
        stop(
        "Package 'chromVAR' is required for addGCBias_ChAI. ",
        "Please install 'chromVAR' to proceed."
        )
    }

    obj1 <- chromVAR::addGCBias(obj1, genome = genome)
    if (any(is.na(SummarizedExperiment::rowData(obj1)$bias))) {
        naList <- is.na(SummarizedExperiment::rowData(obj1)$bias)
        
        if (verbose) {
        warning(paste(sum(naList), "NaNs found within GC Bias", sep = " "))
        }
        
        SummarizedExperiment::rowData(obj1)$bias[which(naList)] <- mean(SummarizedExperiment::rowData(obj1)$bias, na.rm = TRUE)
    }
    return(obj1)
}

#' @title Reformats raw output from makeChromVAR
#'
#' @description \code{reformatChromVAR} Takes a celltype list of chromVARDeviations or a single chromVARDeviations object over all cell types and reformats it into a ChAI object where each assay is a cell type
#'
#' @param chromVARList The output of makeChromVAR, either as a cell type list (withinCellType and exportList are TRUE, which produces a a named list of ChromVAR objects for each cell type, or a single ChromVAR object run across multiple cell types when withinCellType is FALSE
#' @param selectDev A boolean, with default to FALSE, that decides whether or not deviations or Z score values will be used in the output object. Default of FALSE chooses Z-scores.
#'
#' @return A SummarizedExperiment of either z-scores or deviations across all cell types, (each assay is a cell type), similar to the format of a MOCHA object for scATAC data, or the output of makePseudobulkRNA. 
#'              This format is necessary for using ChromVAR with modeling and association functions (model_General, scATAC_Associations, scRNA_Associations, general_Associations)
#'
#' @keywords data_import
#'
#' @export

reformatChromVAR <- function(chromVARList, selectDev = FALSE){

    CellType <- NULL

    if(!methods::is(chromVARList, 'list') & !methods::is(chromVARList, 'chromVARDeviations')){
        stop('chromVARList is not a list or an individual chromVARDeviation.')
    }else if(methods::is(chromVARList, 'chromVARDeviations')){

        allMeta <-  SummarizedExperiment::colData(chromVARList)
        
        tmpList <- reformatMetaData(allMeta)
        
        sampleData <- tmpList[[1]][sort(rownames(tmpList[[1]])),]
        summarizedData <- tmpList[[2]] 
        
        if(selectDev){
            assayType = 'deviations'
        }else{
            assayType = 'z'
        }
        
        fullMat <- SummarizedExperiment::assays(chromVARList)[[assayType]]
        rownames(summarizedData) <- gsub(" ","_",rownames(summarizedData))
        cellNames <- gsub("__.*","", colnames(fullMat))
       
        assayList <- lapply(unique(cellNames), function(XX){
            
                tmpMat <- fullMat[,cellNames == XX]
                newCellType <- paste(gsub(" |//.","_", XX), collapse ="|")
                colnames(tmpMat) <- sub("__","",gsub(paste("^",newCellType, sep=''),"", colnames(tmpMat)))
                tmpMat[,sort(colnames(tmpMat))]
            
            })
        names(assayList) <- unique(cellNames)
        
    }else if(!all(unlist(lapply(chromVARList, class)) == 'chromVARDeviations')){
        
        stop('Some or all indices of chromVARList are not chromVarDeviations objects. Check your list.')
        
    }else{

        allMeta <- do.call('rbind', lapply(chromVARList, SummarizedExperiment::colData))
        #Identify the end of the original metadata, after which CellType, and various other CellType-Sample metadata was tacked on. 
       
        tmpList <- reformatMetaData(allMeta)

        sampleData <- tmpList[[1]][sort(rownames(tmpList[[1]])),]
        summarizedData <- tmpList[[2]] 

        # Pull out the assays needed and rename them after the celltypes. 
        if(selectDev){
            assayType = 'deviations'
        }else{
            assayType = 'z'
        }
        
        assayList <- lapply(seq_along(chromVARList), function(x){
            tmpMat <- SummarizedExperiment::assays(chromVARList[[x]])[[assayType]]
            newCellType <- paste(gsub(" |//.","_", names(chromVARList)[x]), collapse ="|")
            colnames(tmpMat) <- sub("__","",gsub(paste("^",newCellType, sep=''),"", colnames(tmpMat)))
            tmpMat
            })
        
    }
        
    outputSE <- SummarizedExperiment::SummarizedExperiment(assayList, colData = sampleData, 
                        metadata = list('summarizedData' = summarizedData,
                                        'History' = paste("reformatChromVARList", utils::packageVersion("ChAI"))))

    return(outputSE)
}


#' @title \code{reformatMetaData}
#'
#' @description \code{allMeta} Takes metadata from a combinedSampleTileMatrix, combined ChromVAR object, or combinedPseudobulk object and splits it back out into a sample-specific metadata data.frame and a SummarizedExperiment containing population-specific metrics (celltype by sample). 
#'
#' @param allMeta A data.frame with a CellType and Sample column, containing metadata across all samples and cell type. From the first column untill the Celltype column, sample-specific metadata is encoded, and from the Celltype column till the end, Population-level metadata is encoded.
#'
#' @noRd

reformatMetaData <- function(allMeta){
    
    . <- newCellType <- CellType <- NULL
    
    limitation1 <-which(colnames(allMeta) == 'CellType')
    #Transform metadata for cell type, and find the sample-level data
    newCellType <- paste(gsub(" |//.","_", unique(allMeta$CellType)), collapse ="|")
    newSample <- sub("__","",sub(paste("^",newCellType,sep=''),"", allMeta$Sample))
    originMeta <- allMeta[, 1:(limitation1-1)]
    originMeta$Sample = newSample
    sampleData <- dplyr::distinct(as.data.frame(originMeta))
    rownames(sampleData) <- sampleData$Sample


    ## Transform the CellType-Sample metadata
    sumData <- as.data.frame(allMeta[,(limitation1+1):length(colnames(allMeta))]) 
    sumDatalist <- lapply(colnames(sumData), function(x){
            tmpMat <- sumData[,x, drop = FALSE]
            tmpMat$CellType = allMeta$CellType
            tmpMat$Sample = newSample
            newTmp <- as.data.frame(tidyr::pivot_wider(tmpMat,id_cols = 'CellType', names_from = 'Sample',
                            values_from = {{x}}))
            rownames(newTmp) <- newTmp$CellType
            dplyr::select(newTmp, !CellType)
    })

    names(sumDatalist) <-   colnames(sumData)     

    summarizedData = SummarizedExperiment::SummarizedExperiment(sumDatalist, colData = sampleData)
    return(list(sampleData, summarizedData))

}

