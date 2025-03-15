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
                                normalize = FALSE) {

    if(any(!c(cellTypeColumn, sampleColumn) %in% colnames(SO@meta.data))){

        stop('Sample or cell type columns are missing from Seurat object. Please verify that the provided column names are correct.')

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
    cl <- parallel::makeCluster(numCores)
    
    subSample_list <- pbapply::pblapply(cl = cl, X= cellTypeList, .subSampleList, 
                                   sampleList1 = sampleList,
                                   sampleColumn1 = sampleColumn, metadata1 = metadata)
    counts_ls <- lapply(subSample_list, function(XX){    
    
                    bothLists <- lapply(XX, function(YY){
                        
                         if(length(YY) > 0){
                                list(unlist(rowSums(as.matrix(counts[,YY]))), 
                                     unlist(rowSums(as.matrix(counts[,YY] > 0))/length(YY)))
                            }else{
                                emptySample
                            }
                        
                        })
    
                  counts <- do.call('cbind', lapply(bothLists, function(z) z[[1]]))
                  percent <- do.call('cbind', lapply(bothLists, function(z) z[[2]]))
        
                colnames(counts) <- colnames(percent) <- sampleList                                     
                list(counts, percent)
        })
                                     
    countList = lapply(counts_ls, function(x) x[[1]])
    percentList = lapply(counts_ls, function(x) x[[2]])  
                       
    names(countList) <- names(percentList) <- cellTypeList
                         
    parallel::stopCluster(cl)
                         
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
                colData = sampleData
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
        
        additionalMetaData <- lapply(additionalCellData, function(x) {
           

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
                "CellCounts" = allCellCounts
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
        txList$GeneSymbol <- suppressWarnings(AnnotationDbi::mapIds(OrgDb, as.character(txList$gene_id), Seurat_format, "ENTREZID"))
        txList <- plyranges::filter(txList, !is.na(GeneSymbol) & GeneSymbol %in% rownames(rnaSE))
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

        #Subset down the read count matrix to just the transcripts we ahve in this database.
        newSE <- rnaSE[rownames(rnaSE)  %in% txList$GeneSymbol, ]
        
        # Repackage the gene-sample matrices, the sampleData, transcript Genomic Ranges, and associated metadata (countInfo) into one SummarizedExperiment object. 
        attachedList <- plyranges::ungroup(txList)
        SummarizedExperiment::rowRanges(newSE) = attachedList[match(rownames(newSE), names(attachedList))]
       
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
                old_mat <- SummarizedExperiment::assays(rnaSE)[[x]]
                old_mat[is.na(old_mat)] = 0
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