#' @title Generate a ChAI scRNA object (Pseudobulk scRNA)
#' 
#' @description \code{makePseuduoulkRNA} pseodubulks a Seurat object by sample and cell type into a SummarizedExperiment, similar to MOCHA

#'
#' @param SO Seurat Object 
#' @param cellTypeColumn The column of the Seurat object with cell type information
#' @param sampleColumn The column of the Seurat object with sample information
#' @param cellPopulations A list of cell type names to extract or the word 'All' (which is default). Determines which cell populations you want to aggregate and analyze from the Seurat object.
#' @param numCores The number of cores to parallelize this operation over. 
#' @param normalize A boolean to determine whether or not to return DESeq2-normalized counts, or raw aggregated counts. Default is TRUE.
#' @param Seurat_format A string, describing the gene name format of the Seurat object. This is used to annotate gene loci, and convert IDs. See documentation for AnnotationDbi::mapIds to identify formats. Default is 'SYMBOL', but 'ENSEMBL' is also common.
#' @param TxDb A Transcript database to be used for identifying gene locations. Must be installed already. Both TxDb and OrgDb must be provided to annote the gene locations.
#' @param OrgDb An Organism database to match up gene names and locations. Must be installed already. Both TxDb and OrgDb must be provided to annote the gene locations.
#' @param verbose Boolean. Default is TRUE, and will print messages. 
#' @return A SummarizedExperiment carrying pseudobulked average expression per 1000 cells for each cell type. 
#'
#' @keywords data_import
#'
#' @export


makePseudobulkRNA <- function(SO, cellTypeColumn, sampleColumn = "Sample", 
                                 cellPopulations = 'All',
                                numCores = 2,
                                 Seurat_format = 'SYMBOL',
                                 TxDb = NULL, 
                                 OrgDb = NULL,
                                normalize = FALSE,
                             verbose = TRUE) {

    if(any(!c(cellTypeColumn, sampleColumn) %in% colnames(SO@meta.data))){

        stop('Sample or cell type columns are missing from Seurat object. Please verify that the provided column names are correct.')

    }
    
    if(numCores > 3){
    
        warning('User requested to multithread over more than 3 cores. More multithreading will lead to greater memory usage')
    }
    
    if(is.null(TxDb) & !is.null(OrgDb)){

        stop('An OrgDb was provided, but no TxDb was provided. Please provide both, or neither')

    }else if(!is.null(TxDb) & is.null(OrgDb)){

        stop('An TxDb was provided, but no OrgDb was provided. Please provide both, or neither')

    }else if(!methods::is(OrgDb, 'OrgDb') & !is.null(OrgDb)){

        stop('The OrgDb provided is not a Organism Database. Please provide an Organism Database or NULL.')

    }else if(!methods::is(TxDb, 'TxDb') & !is.null(TxDb)){

        stop('The TxDb provided is not a Transcript Database. Please provide a Transcript Database or NULL.')
    }
        
     
    #Generate sample and cell type column in metadata for pseudobulking (after filtering down)
    metadata <- as.data.frame(SO@meta.data)
    counts <- SO@assays$RNA@counts
    if(any(dim(counts)==0)){

        counts <- SO@assays$RNA@counts
    }
    if(any(dim(counts)==0)){

        stop('Error: Cannot find count matrix within Seurat Object')
    }
    metadata$CellTypeColumn = factor(unlist(metadata[,cellTypeColumn]))
    
    if(any(is.na(SO@meta.data[[sampleColumn]]))){
        stop('Some values within the sampleColumn are NA. Please correct the NAs or choose a different column for sample identification.')
    }
    if(all(rownames(metadata) == c(1:dim(metadata)[1]))){

        stop('Seurat metadata has corrupted rownames. Please make sure the Seurat object has cell names saved.')
        
    }

    if(all(tolower(cellPopulations) == 'all')){
        cellPopulations <- unique(metadata$CellTypeColumn)
    }else if(all(cellPopulations %in% unique(metadata[,cellTypeColumn]))){
        metadata <- dplyr::filter(metadata, !!as.name(cellTypeColumn) %in% cellPopulations)
        counts <- counts[,rownames(metadata)]
    }else{
        stop('Within the cellTypeColumn, not all cellPopulations could be found.')
    }

    cellTypeList <- as.character(unique(metadata$CellTypeColumn))
    sampleList <- unique(SO@meta.data[,sampleColumn])
    emptySample = counts[,1]
    emptySample[TRUE] = 0

    if(verbose){     message("Generating sample List")    }
    
    metadata$Cells = rownames(metadata)
    cellTypeListDFs <- dplyr::group_split(metadata, CellTypeColumn)
    names(cellTypeListDFs) = unlist(lapply(cellTypeListDFs, function(XX) unique(XX$CellTypeColumn)))
    ## Subset down to desired populations
     cellTypeListDFs =  cellTypeListDFs[cellPopulations]
    # Create list to iterate over. 
    subSampleList <- pbapply::pblapply(cl = NULL, X =  cellTypeListDFs, function(XX){
                            subMeta = dplyr::group_split(XX, !!as.name(sampleColumn))
                            subCells = lapply(subMeta, function(ZZ) ZZ$Cells)
                            names(subCells) = unlist(lapply(subMeta, function(ZZ) unique(ZZ[,sampleColumn])))
                            subCells
                        })
    
    if(verbose){  message("Processing total counts and percent detection.")    }
                                                     
    if(numCores > 1){
        cl <- parallel::makeCluster(numCores)
    }else{
        cl = NULL
    }
                                                     
    iterList = lapply(subSampleList, function(ZZ){
                    list('Cells' = ZZ, 
                         'Counts' = as.matrix(counts[,unlist(ZZ)]),
                         'sampleList' = sampleList)
                })
    gc()
    rm(counts)
    countList <- pbapply::pblapply(cl = cl, X = iterList, 
                                   .processCounts, returnCounts = TRUE,
                                   emptySample = emptySample)
    gc()
    percentList <- pbapply::pblapply(cl = cl, X = iterList, .processCounts, 
                                     returnCounts = FALSE, emptySample = emptySample)
                       
    names(countList) <- names(percentList) <- cellTypeList
                         
    if(numCores > 1){
        parallel::stopCluster(cl)
    }
    rm(iterList)
                         
    ##Sort sample-level metadata
    cellColDataNoNA <- BiocGenerics::Filter(function(x) {
        !all(is.na(x))
    }, SO@meta.data)   

    sampleSpecificColumns <- dplyr::group_by(cellColDataNoNA, !!as.name(sampleColumn))  
    sampleSpecificColumns <- dplyr::summarize(sampleSpecificColumns, dplyr::across(dplyr::everything(), dplyr::n_distinct)) 
    sampleSpecificColumns <- dplyr::select(sampleSpecificColumns, tidyselect::where(~all(.x==1)))
    sampleData <- dplyr::distinct(cellColDataNoNA[,c(sampleColumn,colnames(sampleSpecificColumns))])
    # Set sampleIDs as rownames
    rownames(sampleData) <- sampleData[[sampleColumn]]
    sampleData$Sample = sampleData[[sampleColumn]]

    ## Save percent detected for each gene.
    percentDetected <- SummarizedExperiment::SummarizedExperiment(
               percentList,
                colData = sampleData[sampleList,]
    )

    #Get cell numbers
    cellCounts <- as.data.frame(table(metadata[, sampleColumn], metadata[, 'CellTypeColumn']))
    names(cellCounts) <- c(sampleColumn, "CellPop", "CellCount")
    cellCounts <- tidyr::pivot_wider(
        cellCounts,
        id_cols = "CellPop",
        names_from = sampleColumn,
        values_from = "CellCount"
    )
    allCellCounts <- as.data.frame(cellCounts[, -1])
    rownames(allCellCounts) <- cellCounts$CellPop
    allCellCounts <- allCellCounts[sort(cellPopulations),]

    #Process the meta data and extract the total counts and features for each population, as well as cell counts. 
    cellColDataCopy <- data.frame(metadata)

    cellColDataCopy[] <- lapply(cellColDataCopy, function(x) {
        utils::type.convert(as.character(x), as.is = TRUE)
    })

    # Assume all numeric columns are to be saved as additionalCellData
    isNumericCol <- unlist(lapply(cellColDataCopy, function(x) is.numeric(x)))
    additionalCellData <- colnames(cellColDataCopy)[isNumericCol]

    # Group by Sample (rows) and cellPop (columns)
    summarizedData_df <- dplyr::group_by(cellColDataCopy, CellTypeColumn, !!as.name(sampleColumn))
    if (!is.null(additionalCellData)) {
        
        if(verbose){     message("Processing additional metadata.")    }
        
        additionalMetaData <- pbapply::pblapply(cl = NULL, X = additionalCellData, function(x) {
           

            suppressMessages(
                summarizedData2 <- dplyr::summarize(summarizedData_df, meanValues = mean(!!as.name(x), na.rm = TRUE, .groups = "drop")) 
            )

            summarizedData2 <- tidyr::pivot_wider(summarizedData2,
                    id_cols = CellTypeColumn,
                    names_from = !!as.name(sampleColumn),
                    values_from = meanValues
                )

            summarizedData2 <- as.data.frame(summarizedData2)
            rownames(summarizedData2) <- summarizedData2[['CellTypeColumn']]
            summarizedData2 <- summarizedData2[, -1, drop = FALSE]

            # Filter to specific cellPopulations
            summarizedData2 <- summarizedData2[
                rownames(summarizedData2) %in% cellPopulations, ,
                drop = FALSE
            ]

            summarizedData2[sort(cellPopulations),sort(colnames(summarizedData2))]
        })
        names(additionalMetaData) <- additionalCellData

    }else if (is.null(additionalCellData)) {

        additionalMetaData <- NULL

    }
                      
    remove(cellColDataCopy)

    summarizedData <- SummarizedExperiment::SummarizedExperiment(
            append(
            list(
                "CellCounts" = allCellCounts[rownames(additionalMetaData[[1]]),
                                             colnames(additionalMetaData[[1]])]
            ),
            additionalMetaData
            ),
            colData = sampleData[sort(rownames(sampleData)),]
        )

    rnaSE <-  SummarizedExperiment::SummarizedExperiment(countList, colData = sampleData,
         metadata = list('summarizedData' = summarizedData,'detectionRates' = percentDetected,
            History = paste("makePseudobulkRNA", utils::packageVersion("ChAI"))))

    #Decide if you are interweaving this data with genomic databases/locations. 
                                  
    #Pull in Transcript and Organism databases. 
    if(!is.null(TxDb) & !is.null(OrgDb)){
        
        newSE <- linkToGenome(rnaSE = rnaSE, gene_format = Seurat_format, 
                              TxDb= TxDb, OrgDb = OrgDb)
       
    }else{
        newSE <- rnaSE

    }

    if(normalize){
        newSE <- normalizePseudobulk(newSE, sampleColumn = sampleColumn)
    }
    return(newSE)
}
                                  
                                  
#' @title \code{subSampleList}
#'
#' @description \code{subSampleList} Helper function for makePseudobulkRNA. 

#'
#' @param cellType  the output of makePseudobulkRNA a
#' @param sampleList1 List of sample names
#' @param sampleColumn1 Name of metadata column with sample names
#' @param metadata1 Seurat metadata
#' @return a list of metadata by celltype and sample
#'
#' @noRd

.subSampleList <- function(cellType, sampleList1, sampleColumn1, metadata1){
        
        cellNameList <- lapply(sampleList1, function(z){
            rownames(dplyr::filter(metadata1, 
                            CellTypeColumn == cellType & !!as.name(sampleColumn1) == z))
        })
        return(cellNameList)
}

#' @title \code{processCounts}
#'
#' @description \code{processCounts} Helper function for makePseudobulkRNA. 

#'
#' @param sampleList_index One index from the subSampleList
#' @param returnCounts Boolean. If true, will return count matrix. If false, will return percent detected. 
#' @param emptySamples A matrix representing an empty sample (all genes 0)
#' @return a list of matrices, either counts or percent detected
#'
#' @noRd
                   

.processCounts <- function(sampleList_index, returnCounts = TRUE, emptySamples){
    
    if(returnCounts){     

    matList <- lapply(sampleList_index[[1]], function(YY){
                        
            if(length(YY) > 1){
                    unlist(rowSums(sampleList_index[[2]][,YY]))
            }else if(length(YY) == 1){
                sampleList_index[[2]][,YY]
            }else{
                emptySample
            }

        }) 

    }else{

        matList <- lapply(sampleList_index[[1]], function(YY){

                if(length(YY) > 1){
                        unlist(rowSums(sampleList_index[[2]][,YY] > 0))/length(YY)
                }else if(length(YY) == 1){
                    as.integer(sampleList_index[[2]][,YY] > 0)
                }else{
                    emptySample
                }

            }) 
	
    }

    mat <- do.call('cbind', matList)
    colnames(mat) <- names(sampleList_index[[1]])
    ## Identify and fill in missing samples with NAs. 
    if(any(!sampleList_index[[3]] %in% colnames(mat))){
        
        missingSamples = sampleList_index[[3]][!sampleList_index[[3]] %in% colnames(mat)]
        
        for(newCol in missingSamples){
        
            mat = cbind(mat, newCol = 0)
            colnames(mat)[colnames(mat) == 'newCol'] = newCol
        }
        
    }
    mat = mat[,sampleList_index[[3]]]
    
    return(mat)
}

#' @title \code{linkToGenome}
#'
#' @description \code{linkToGenome} Helper function for makePseudobulkRNA. Links genes to genomic loci, and tosses poorly annotated transcripts. 
#'
#' @param rnaSE the output of makePseudobulkRNA
#' @param gene_format A string, describing the gene name format of the Seurat object. This is used to annotate gene loci, and convert IDs. See documentation for AnnotationDbi::mapIds to identify formats. Default is 'SYMBOL', but 'ENSEMBL' is also common.
#' @param TxDb A Transcript database to be used for identifying gene locations. Must be installed already. Both TxDb and OrgDb must be provided to annote the gene locations.
#' @param OrgDb An Organism database to match up gene names and locations. Must be installed already. Both TxDb and OrgDb must be provided to annote the gene locations.
#' @return a list of metadata by celltype and sample
#'
#' @keywords data_import
#'
#' @export

linkToGenome <- function(rnaSE, TxDb= NULL, OrgDb = NULL, gene_format = 'SYMBOL'){

    ## fix global bindings
    ALIAS <- SYMBOL <- NULL
    if(!is.null(TxDb) & !is.null(OrgDb)){
        if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
                stop(
                "Package 'AnnotationDbi' is required for adding gene location information.",
                "Please install AnnotationDbi'from Bioconductor to proceed."
                )
        }

        #Extact the GenomicRanges for all genes, while filtering out non-standard chromosomes that mess things up. 
        txList <- suppressMessages(suppressWarnings(S4Vectors::stack(GenomicFeatures::genes(TxDb,
                                                single.strand.genes.only = FALSE))))
        txList <- suppressMessages(GenomeInfoDb::keepStandardChromosomes(sort(txList), 
                                    species='Homo_sapiens',
                                        pruning.mode = 'coarse'))
        #Reduce ranges merges transcripts for the same gene into one, so that we can one ball-park stop and end. 
        txList <- plyranges::reduce_ranges_directed(plyranges::group_by(txList, gene_id))
        txList$GeneSymbol <- suppressWarnings(AnnotationDbi::mapIds(OrgDb, as.character(txList$gene_id), gene_format, "ENTREZID"))
        txList <- plyranges::filter(txList, !is.na(GeneSymbol))
        txList$MultipleMappings = txList$GeneSymbol %in% txList$GeneSymbol[duplicated(txList$GeneSymbol)]
        txList <- plyranges::slice(plyranges::group_by(txList, GeneSymbol), 1)
        names(txList) <- txList$GeneSymbol
        txList <- plyranges::ungroup(txList)
        promoters <- GenomicRanges::promoters(txList)
        promoters <- paste(GenomicRanges::seqnames(promoters), ":", GenomicRanges::start(promoters),"-", GenomicRanges::end(promoters), sep = '')
        txList$Promoters = promoters
        txList$Strand_Direction = as.character(GenomicRanges::strand(txList))
        ## Some genes are mapped to both positive and negative strands, and some miRNAs are provisionally mapped to multiple locations. 
        ## Let's just take the first location for these. 
        
        if(any(!rownames(rnaSE)  %in% txList$GeneSymbol)){

            ## if it's a gene symbol, sometimes there's mismatches, with multiple potential gene names. 
            ## let's make sure that isn't the case, and match up names that don't match the database for SYMBOL. 
            if(gene_format == 'SYMBOL'){
               lostTxs <- rownames(rnaSE)[!rownames(rnaSE)  %in% txList$GeneSymbol]
                #Check for missed due to '.2' gene name versions
               
               aliasAnnot <- AnnotationDbi::select(OrgDb, keys=lostTxs,
                        columns=c("SYMBOL","ENTREZID"), 
                                                   keytype="ALIAS")
                # Filter to matching SYMBOL
                # Identify all ALIASes that were matched well
                aliasAnnot1  <- dplyr::filter(aliasAnnot , !is.na(SYMBOL))
                # Identify all ALIASes that didn't. 
                # Remove .[0-9] and try again.
                aliasAnnot2  <- dplyr::filter(aliasAnnot , is.na(SYMBOL))
                aliasAnnot2  <- dplyr::mutate(aliasAnnot2, ALIAS2 =
                                              gsub("\\.[0-9].*","", ALIAS))
                aliasAnnot3 <- AnnotationDbi::select(OrgDb, 
                            keys=unique(aliasAnnot2$ALIAS2),
                        columns=c("SYMBOL","ENTREZID"), 
                                                   keytype="ALIAS") 
                aliasAnnot3 <- dplyr::filter(aliasAnnot3, !is.na(SYMBOL))
                aliasAnnot3 <- dplyr::inner_join(
                                    aliasAnnot2[,c('ALIAS', 'ALIAS2')],
                                    aliasAnnot3, 
                            by = c('ALIAS2' = 'ALIAS'))[,
                                        c('ALIAS','SYMBOL','ENTREZID')]
                aliasAnnot4 <- dplyr::filter(rbind(aliasAnnot1, aliasAnnot3),
                                             SYMBOL  %in% txList$GeneSymbol &
                                            !SYMBOL %in% rownames(rnaSE))
                aliasAnnot5 <- dplyr::slice_head(
                                    dplyr::group_by(aliasAnnot4, ALIAS),
                                    n=1)
              
                subTx <- txList[txList$GeneSymbol %in% aliasAnnot5$SYMBOL]
                subTx$Alias = subTx$GeneSymbol
                subTx$GeneSymbol = aliasAnnot4$ALIAS[match(subTx$GeneSymbol, 
                                                       aliasAnnot5$SYMBOL)]
                names(subTx) = subTx$GeneSymbol
                subTx$Alias = aliasAnnot4$ALIAS[match(subTx$GeneSymbol, 
                                                       aliasAnnot5$SYMBOL)]
                txList = sort(c(txList, subTx))

              
            }
            #Identify Transcripts that were removed.
            lostTxs <- rownames(rnaSE)[!rownames(rnaSE)  %in% txList$GeneSymbol]
            
            #Subset down the read count matrix to just the transcripts we have in this database.
            newSE <- rnaSE[rownames(rnaSE)  %in% txList$GeneSymbol, ]
            
            message(stringr::str_interp("${length(lostTxs)} transcripts were not found within the provided transcript/organism database, or               did not have a clear genomic location within the standard chromosome assembly for your organism. 
                Filtered transcript IDs are saved in metadata, via obj@metadata[['filteredIDs']]."))
                    
            newSE@metadata = append(newSE@metadata, list('filteredIDs' = lostTxs))

        }
            
        #Subset down the read count matrix to just the transcripts we have in this database.
        newSE <- rnaSE[rownames(rnaSE)  %in% txList$GeneSymbol, ]
        
        # Repackage the gene-sample matrices, the sampleData, transcript Genomic Ranges, and associated metadata (countInfo) into one SummarizedExperiment object. 
        attachedList <- plyranges::ungroup(txList)
        SummarizedExperiment::rowRanges(newSE) = attachedList[match(rownames(newSE), names(attachedList))]
       
    }else{
        newSE <- rnaSE
        warning('No TxDb and/or OrgDb provided.')
    }
    return(newSE)
}



#' @title Normalize a ChAI scRNA object using DESEq2
#'
#' @description \code{normalizePseudobulk} Takes the output of makePseudobulkRNA and normalizes it. 

#'
#' @param rnaSE  the output of makePseudobulkRNA a
#' @param sampleColumn The column of the Seurat object with sample information
#' @return A SummarizedExperiment with normalized average expression
#'
#' @keywords data_import
#'
#' @export

normalizePseudobulk <- function(rnaSE, sampleColumn = 'Sample'){

    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop(
        "Package 'DESeq2' is required for normalizePseudobulk. ",
        "Please install 'DESeq2' to proceed."
        )
    }

    allMat <- lapply(names(SummarizedExperiment::assays(rnaSE)), function(x){
                print(x)
                old_mat <- SummarizedExperiment::assays(rnaSE)[[x]]
                old_mat[is.na(old_mat)] = 0
                old_mat = round(old_mat) #ensure integers
                suppressMessages(
                    dds <- DESeq2::DESeqDataSetFromMatrix(countData = old_mat[,colSums(old_mat) > 0],
                              colData = SummarizedExperiment::colData(rnaSE)[colSums(old_mat) > 0,],
                                design = stats::as.formula(paste( '~', sampleColumn)))
                )

                suppressMessages(dds <- DESeq2::estimateSizeFactors(dds))
                suppressMessages(new_mat <- DESeq2::counts(dds, normalize = TRUE))
                if(any(colSums(old_mat) == 0)){
                    filled_data <- do.call('cbind', lapply(which(colSums(old_mat) == 0), function(x){
                            rep(0, dim(new_mat)[1])
                    }))
                    new_mat <- cbind(new_mat, filled_data)
                    
                }
                
                new_mat <- new_mat[,sort(colnames(new_mat))]
        })

    names(allMat) <- names(SummarizedExperiment::assays(rnaSE))
    
    se <- SummarizedExperiment::SummarizedExperiment(allMat, 
                        colData =SummarizedExperiment::colData(rnaSE)[sort(colnames(rnaSE)),],
                        metadata = rnaSE@metadata)
    SummarizedExperiment::rowRanges(se) <-  SummarizedExperiment::rowRanges(rnaSE)

    se@metadata$History <- append(se@metadata$History, paste("normalizePseudobulk", utils::packageVersion("ChAI")))
    return(se)
    
    }

CellTypeColumn <- gene_id <- GeneSymbol <- meanValues <- NULL