#' @title Combine multiple Association Models into a network of interaction features.
#'
#' @description \code{getRegulation} takes multiple associative models and filters down to links
#'                That show a consistent relationship that is significant across associations. 
#'                In other words, for a given TF to bind a tile and influence gene expression, the motif should be significantly
#'                  associated with gene expression and tile accessibility, and the gene expression should be significantly related to tile accessibility.\
#'                           
#'
#'
#' @param AssociationModelList A list of ChAI association models. Must include at least 3 models. 
#' @param modelAnnotations An optional list of strings describing each association model. Each string would be in the form Modality1_Modality2, and match the order of the association used. This is useful if you ran multiple general modalities. For example, modelAnnotations could be list("ATAC_RNA", "RNA_Protein", "RNA_ChromVAR", "ATAC_Protein", "ATAC_ChromVAR") for an AssociationModelList of list(GeneTile_Association, RNA_Association, RNA_Association, scATAC_Association, scATAC_Association). The order of the modalities must reflect the order of the association model - i.e. scATAC is always modeled using other modalities as a variable, and scRNA is always modeled using other modalities as a variable in the function (except for scATAC). If no modelAnnotations is given, it will simply number the association models 1,2,3,4, etc...
#' @param interactionDirection An optional list of the interaction direction, whicb can be useful if you want to export tables for import into cytoscape or other network visualization softwares. Using this list, it will label one modality as a source and another as a target. For example, we expect a tile's accessibility to potentially regulate a gene's expression, a transcription factor to regulate a tile accessibility, and an olink protein's concentration to influence a transcription factor's activity. So for an AssociationModelList with the model annotations of list('ATAC_RNA', 'ATAC_ChromVAR', 'ChromVAR_Olink'), we would use the interactionDirection list of list('Source_Target', 'Target_Source', 'Target_Source') to describe the direction of the regulation we expect.
#' @param FDR_threshold A number between 0 and 1, used as a significance threshold for associations. Any pairwise associations between modalities that have an FDR below this threshold are labeled as significant and used to build the network. 
#' @param annotation_threshold Default is NULL, in which case it will be set to the length(AssociationModelList)-1, which will require 
#'              All the modalities to be significantly associated with each other (i.e. A TF is sig associated to a tile and a gene, and those tiles and genes are also significantly associated with each other). If you want a less stringent threshold, you can set it to less than that manually. 
#' @param controlFactors A optional list of strings that describe various factors we want to remove as sources of noise. These factors should ideally represent a strong technical artifact of somekind (like sequencing depth or other). If a given association model for two features is significant, it will be removed from downstream analysis. For example, if there is a significant interaction between IFNG protein levels and IRF1 ChromVAR activity, that link should be kept in the network. But if IFNG protein level is also strongly associated with limit of detection for across batches, then removing that link might make sense. You don't know if it's a significant interaction or a technical noise. This decision should be made after looking at each individual model association as well. For associations using scATAC, ZI_FragNumber should be used because drop-out will impact the degree of zero-inflation, so any significant Zero-inflated interaction term (like ZI_Age), but not continuous term, may be driven by technical noise and not the variable of interest, like Age. 
#' @param backgroundThreshold A number between 0 and 1 that is used as the FDR threshold to determine whether a controlFactor is significant, and thus a given pairwise association should be removed. 
#' @param cytoscapeOutput A boolean. Default is FALSE, in which case it will merge significant links together into one data.frame such that all features that pass the annotation_threshold are included. If set to TRUE, it will return a list of two data.frames, which are node and edge tables that can be imported into Cytoscape or other network visualization softwares. 
#' @return a dataframe or dataframes of potential networks
#'
#' @importFrom data.table :=
#'
#' @export
#' @keywords networks
getRegulation <- function(AssociationModelList, 
                          modelAnnotations = NULL,
                          interactionDirection = NULL,
                          FDR_threshold = 0.1,
                          annotation_threshold= NULL,
                          controlFactors = 'ZI_FragNumber',
                          backgroundThreshold=  0.1, 
                          cytoscapeOutput = FALSE) {
    
    ifelse(length(AssociationModelList) <= 2, stop('AssociationModelList must include three or more models. You have provieded fewer than three models.'), "pass")
    if(is.null(annotation_threshold)){
        annotation_threshold = length(AssociationModelList) - 1
    }else if(annotation_threshold > length(AssociationModelList)-1 | annotation_threshold < 2){
        
           stop('The annotation_threshold represents the number of times a given feature must be significantly associated with another modality. For example, a Gene must be significantly related to both a TF activity and a Tile accessibility (2 modalities) to be kept in a Gene Regulatory network. This number can be as low as 2, or as high as n-1, where n is the length of AssociationModelList. You have provided an annotation_threshold that is below 2 or greater than n-1, which will not work.')
        
    }

    #Verify object types       
    if(any(unlist(lapply(AssociationModelList, function(x) !methods::is(x,'SummarizedExperiment'))))){
        
        stop('At least one of association models in the AssociationModelList is not a SummarizedExperiment.')
    }
                         
    if(any(unlist(lapply(AssociationModelList, function(x) !any(names(x@metadata) == 'Type'))))){
        
        stop('At least one of the inputs (GeneMotifLinks, GeneTileLinks, or TileMotifLinks) are not SummarizedExperiment generated by ChAI.')
    }
                         

    modelTypes <-   unlist(lapply(AssociationModelList, function(x) x@metadata$Type))
                                  
    lapply(AssociationModelList, isChAIObject, type = 'model')
                                  
                                                  
    if(any(duplicated(modelTypes)) & !is.null(modelAnnotations)){

        UniqueAssociations = FALSE
        modelTypeList = c(gsub("_.*", "", unlist(modelAnnotations)), gsub(".*_", "", unlist(modelAnnotations)))
        #Each type of model should be duplicated at 
        if(!all(as.data.frame(table(modelTypeList))$Freq >= 2)){
        
            stop("All modalities must be associated with at least two other modalities to form a network.")
            
        }

        
    }else if(any(duplicated(modelTypes) & is.null(modelAnnotations))){
        
        UniqueAssociations = FALSE
        modelAnnotations = gsub("^sc", "", gsub('_Associations', '_Other', modelTypes))
        
        if(!cytoscapeOutput){
            dupTypes <- unique(modelTypes[ which(duplicated(modelTypes))])
            stop(stringr::str_interp("You have provided an AssociationModelList that includes multiple association models of the same type, and you have not provided any modelAnnotations. Under these conditions, getRegulation does not know which modalities to bind. Duplicated model types are ${paste(dupTypes, collapse=', ')}."))
            
        }
        
        
    }else if(!is.null(modelAnnotations)){
        
        UniqueAssociations = TRUE
        
        modelTypeList = c(gsub("_.*", "", unlist(modelAnnotations)), gsub(".*_", "", unlist(modelAnnotations)))
        #Each type of model should be duplicated at 
        if(!all(as.data.frame(table(modelTypeList))$Freq >= annotation_threshold)){
        
            stop("All modalities must be associated with at least ${annotation_threshold} other(s) to form a network.")
            
        }
        
    }else{
        
        UniqueAssociations = TRUE
        modelAnnotations = gsub("^sc", "", gsub('_Associations', '_Other', modelTypes))
    }

        
    ## If you want a cytoscape output, it will annotate the source and target of the interaction you are testing
    ## from the modelAssociationList, if you give it an interactionDirection. 
    #If you don't, it'll just label it Modality1 and Modality2, instead of source and taret. 
    if(cytoscapeOutput & !is.null(interactionDirection)){
    
        if(!all(interactionDirection %in% c('source_target','target_source'))){
                stop(stringr::str_interp("interactionDirection must be list that only contains 'source_target' or 'target_source'"))
        }
        
    }else if(cytoscapeOutput){
        
        interactionDirection = rep(c('Modality1_Modality2'), length(modelAnnotations))
    }
                  
                                  
    #Verify that controlFactors are present in at least one model
    names(AssociationModelList) <- modelTypes
        
    assayNames <- lapply(AssociationModelList, function(x) { names(SummarizedExperiment::assays(x)) })                          
    if(!is.null(controlFactors)){
        
        if(!all(controlFactors %in% unlist(assayNames))){
            controlFactorMissing <- controlFactors[which(!controlFactors %in% unlist(assayNames))]
            stop(stringr::str_interp("These controlFactors were not found within any models: ${paste(controlFactorMissing, collapse=', ')}"))

        }
        
    }else{
        
        controlFactors = character()
        
    }
                
    ## Get estimates from each object. If one object doesn't have any significant assocations, the throw an error. 
    AllEstimates <- lapply(seq_along(modelTypes), function(x){
        
           
            associationModel <- AssociationModelList[[x]]
        
            leftOverFactors = controlFactors[controlFactors %in% assayNames[[x]]]
        
            ## If a factor has a ZI_ term, it will only removed from a ZI_ portion. 
            ## This function will also pull out ZI_exp2 for all zero-inflated portions. 
            
            #Identify non-zero inflated associations   
            specControls = grep("ZI_", controlFactors[controlFactors %in% assayNames], value = TRUE, invert= TRUE)
            if(length(specControls) == 0){ specControls <- NULL}
        
            interactionFactors = grep("ZI_", grep('exp2', names(SummarizedExperiment::assays(AssociationModelList[[x]])), value = TRUE), 
                                                    value = TRUE, invert = TRUE)
            if(length(interactionFactors) == 0 ){
                    stop(stringr::str_interp('No interaction term (exp2) was found for index ${x} of the AssociationModelList.'))
            }
                         
            baseEstimates <- do.call('rbind', lapply(interactionFactors, function(ZZ){
                    baseEstimatesTmp <- getEstimates(modelSE = AssociationModelList[[x]], factor = ZZ, FDR_threshold = FDR_threshold,
                                                      filterFactors = specControls, backgroundThreshold = backgroundThreshold)
                    if(any(dim(baseEstimatesTmp) ==0)){
                           stop(stringr::str_interp('Association Model number ${x} has no significant interactions'))
                    }
                    baseEstimatesTmp$Class = 'Continuous'
                    baseEstimatesTmp
                }))
        
            if(any(grep('ZI_exp2',assayNames[[x]]))){

                #Identify zero-inflated terms
                specControlsZ = grep("ZI_", controlFactors[controlFactors %in% assayNames], value = TRUE)
                if(length(specControlsZ) == 0){ specControlsZ <- NULL}
                
                interactionFactors_ZI = grep('ZI_exp2', names(SummarizedExperiment::assays(AssociationModelList[[x]])), value = TRUE)
                
                ZI_Estimates <- do.call('rbind', lapply(interactionFactors_ZI, function(ZZ){
                    
                    ZI_EstimatesTmp <- getEstimates(modelSE = AssociationModelList[[x]], factor = ZZ, FDR_threshold = FDR_threshold,
                                                    filterFactors = specControlsZ, backgroundThreshold = backgroundThreshold)
                    if(dim(ZI_EstimatesTmp)[1] != 0){
                        
                        ZI_EstimatesTmp$Class = 'ZeroInflated'
                        ZI_EstimatesTmp
                    }else{
                    
                        NULL
                    }
                    
                }))
                

                return(rbind(baseEstimates, ZI_Estimates))

            }else{
                return(baseEstimates)
            }
        
        })
    
    ## Make sure all the models have significant associations.
    if(any(unlist(lapply(AllEstimates, function(x) dim(x)[1] == 0)))){
        failNum <- which(unlist(lapply(AllEstimates, function(x) dim(x)[1] == 0)))
        stop(stringr::str_interp('Association Model number ${failNum} has no significant interactions'))
    }
                                       
    ## Now we need to filter down to features that showed up in at least n-1 association models, 
                                       #which means they were significantly associated with every other modality tested. 
          
    #Split the model annotation list into a list of vectors, where the first index of each vector is the first modality, and the second index is the second.
    annotationList <- lapply(modelAnnotations, function(XX){
        
            unlist(strsplit(x= ifelse(grepl('GeneTile', XX), c('ATAC_RNA'), XX), split = "_"))
        
        })
                                       
    featureList <- do.call('rbind', lapply(seq_along(AllEstimates), function(XX) {  

        if(grepl('GeneTile', modelTypes[[XX]])){
                          
            GeneList <- unique(AllEstimates[[XX]]$Genes)
            TileList <- unique(AllEstimates[[XX]]$Tiles)
            data.frame(Feature = c(TileList, GeneList), 
                       Type = c(rep(annotationList[[XX]][1], length(TileList)), rep(annotationList[[XX]][2], length(GeneList)))) 
             
        }else if(grepl('scATAC_', modelTypes[[XX]])){
            
            GeneralList <- unique(AllEstimates[[XX]]$General)
            TileList <- unique(AllEstimates[[XX]]$Tiles)
            data.frame(Feature = c(TileList, GeneralList), 
                       Type = c(rep(annotationList[[XX]][1], length(TileList)), rep(annotationList[[XX]][2], length(GeneralList)))) 
             
        }else if(grepl('scRNA_', modelTypes[[XX]])){
            
            GeneList <- unique(AllEstimates[[XX]]$Genes)
            GeneralList <- unique(AllEstimates[[XX]]$General)
            data.frame(Feature = c(GeneList, GeneralList), 
                       Type = c(rep(annotationList[[XX]][1], length(GeneList)), rep(annotationList[[XX]][2], length(GeneralList)))) 
        
        }else if(grepl('ChromVAR_', modelTypes[[XX]])){
            
            TFList <- unique(AllEstimates[[XX]]$TFs)
            GeneralList <- unique(AllEstimates[[XX]]$General)
            data.frame(Feature = c(TFList, GeneralList), 
                       Type = c(rep(annotationList[[XX]][1], length(TFList)), rep(annotationList[[XX]][2], length(GeneralList)))) 
        
        }else if(grepl('General', modelTypes[[XX]])){
            
            GeneralList1 <- unique(AllEstimates[[XX]]$Modality1)
            GeneralList2 <- unique(AllEstimates[[XX]]$Modality2)
            data.frame(Feature = c(GeneralList1, GeneralList2), 
                       Type = c(rep(annotationList[[XX]][1], length(GeneralList1)), rep(annotationList[[XX]][2], length(GeneralList2)))) 
        
        } else {
             
             stop("Model type not recognized. Model type should have already been checked, so some internal error likely happened.")
        
        }
    
    }))
            
                               
    filtFeatureList <- dplyr::summarize(dplyr::group_by(featureList, Feature, Type), Number = dplyr::n())
    
    subFeatureList <- dplyr::filter(filtFeatureList, Number >= annotation_threshold)
        
    if(any(dim(subFeatureList) == 0)){
    
        stop("Features are not consistently association across modalities. Ensure you have provided every available combination of data types, or lower the annotation_threshold set too high.")
        
    }
    
    ## We need to figure out how to bind these matrices together, while keeping the Estimate value and FDR intact for latter manipulation. 
    ## Ideally, there's only one model object for each type of associations (GeneTile, scATAC, scRNA, chromVAR, and General), but this isn't always necessarily true.
    ## Especially if you want to include scRNA, scATAC, Olink, and ChromVAR. You would then have two scATAC_General - one for scATAC-Olink and one for scATAC-ChromVAR
    ## If there is, then we can simply add a flag in front - Estimate becomes GeneTile_Estimate or scATAC_General_Estimate, for example. 
    ## If there isn't, then we can simply add a number for the index of the list, like Estimate1 for the first. 
                                       
    ## We may not always want to run with Olink because some TFs repond to steroids which won't appear in Olink. 
                                       
    ### Now let's join each association estimate data.frame together. 
                                                            
    keptEstimates <- lapply(seq_along(AllEstimates), function(XX){
        
        tmpMat <- AllEstimates[[XX]]
        
        if(grepl('GeneTile', modelTypes[[XX]])){
                            
            tmpEstimates <- dplyr::select(tmpMat, Estimate, FDR, Tiles, Genes, Class)
            TileList = dplyr::filter(subFeatureList, Type == annotationList[[XX]][1])
            GeneList = dplyr::filter(subFeatureList, Type == annotationList[[XX]][2])
            tmpEstimates <- dplyr::filter(tmpEstimates, Tiles %in% TileList$Feature, Genes %in% GeneList$Feature)
             
            Modality1 = 'Tiles'
            Modality2 = 'Genes'
             
        }else if(grepl('scATAC_', modelTypes[[XX]])){
            
            tmpEstimates <- dplyr::select(tmpMat, Estimate, FDR, Tiles, General, Class)
            TileList = dplyr::filter(subFeatureList, Type == annotationList[[XX]][1])
            GeneralList = dplyr::filter(subFeatureList, Type == annotationList[[XX]][2])
            tmpEstimates <- dplyr::filter(tmpEstimates, Tiles %in% TileList$Feature, General %in% GeneralList$Feature)
             
            Modality1 = 'Tiles'
            Modality2 = 'General'
             
        }else if(grepl('scRNA_', modelTypes[[XX]])){
            
            tmpEstimates <- dplyr::select(tmpMat, Estimate, FDR, Genes, General, Class)
            GeneList = dplyr::filter(subFeatureList, Type == annotationList[[XX]][1])
            GeneralList = dplyr::filter(subFeatureList, Type == annotationList[[XX]][2])
            tmpEstimates <- dplyr::filter(tmpEstimates, Genes %in% GeneList$Feature, General %in% GeneralList$Feature)
             
            Modality1 = 'Genes'
            Modality2 = 'General'
        
        }else if(grepl('ChromVAR_', modelTypes[[XX]])){
            
            tmpEstimates <- dplyr::select(tmpMat, Estimate, FDR, General, TFs, Class)
            TFList = dplyr::filter(subFeatureList, Type == annotationList[[XX]][1])
            GeneralList = dplyr::filter(subFeatureList, Type == annotationList[[XX]][2])
            tmpEstimates <- dplyr::filter(tmpEstimates, TFs %in% TFList$Feature, General %in% GeneralList$Feature)
             
            Modality1 = 'TFs'
            Modality2 = 'General'
        
        }else if(grepl('General', modelTypes[[XX]])){
            
            tmpEstimates <- dplyr::select(tmpMat, Estimate, FDR, Modality1, Modality2, Class)
            GeneralList1 = dplyr::filter(subFeatureList, Type == annotationList[[XX]][1])
            GeneralList2 = dplyr::filter(subFeatureList, Type == annotationList[[XX]][2])
            tmpEstimates <- dplyr::filter(tmpEstimates, Modality1 %in% GeneralList1$Feature, Modality2 %in% GeneralList2$Feature)
             
            Modality1 = 'Modality1'
            Modality2 = 'Modality2'
        
        } else {
             
             stop("Model type not recognized. Model type should have already been checked, so some internal error likely happened.")
        
        }
             
        if(dim(tmpEstimates)[1] == 0){
        
            stop(stringr::str_interp("Index ${XX} of the AssociationModelList doesn't have any significant features that overlapped with other association models. Please double check your AssociationModelList and/or the modelAnnotations you provided. If you mislabeled or inconsistently labeled your modelAnnotations, then the function won't identify if a Feature is significant associated with multiple modalities. The order of the label matters. The modality being modeled should come first, and the modality acting as a factor in the model should come second. For example 'scRNA_Other' or 'scATAC_scRNA'. The label that fails here is ${modelAnnotations[[XX]]}."))
            
        }
           
        if(UniqueAssociations & !cytoscapeOutput){
            
            colnames(tmpEstimates)= gsub(Modality1, annotationList[[XX]][1], colnames(tmpEstimates))
            colnames(tmpEstimates)= gsub(Modality2, annotationList[[XX]][2], colnames(tmpEstimates))
            
            colnames(tmpEstimates)[grepl('Estimate|FDR|Class', colnames(tmpEstimates))] = paste(modelAnnotations[XX], 
                                                                                    grep('Estimate|FDR|Class', colnames(tmpEstimates), value = TRUE),
                                                                                          sep ='_')
            
        }else if(!cytoscapeOutput){
            
            colnames(tmpEstimates)= gsub(Modality1, annotationList[[XX]][1], colnames(tmpEstimates))
            colnames(tmpEstimates)= gsub(Modality2, annotationList[[XX]][2], colnames(tmpEstimates))
            
            colnames(tmpEstimates)[grepl('Estimate|FDR|Class', colnames(tmpEstimates))] = paste(grep('Estimate|FDR|Class', 
                                                                                            colnames(tmpEstimates), value = TRUE), XX, sep ='_')
            
        }else{

            #Standardize the format for cytoscape outputs (Edge table and node table)
            Modality1_Direction = gsub("_.*","", interactionDirection[[XX]])
            Modality2_Direction = gsub(".*_","", interactionDirection[[XX]])
            
            colnames(tmpEstimates)= gsub(Modality1, Modality1_Direction, colnames(tmpEstimates))
            colnames(tmpEstimates)= gsub(Modality2, Modality2_Direction, colnames(tmpEstimates))
            
            tmpEstimates[paste(Modality1_Direction,'_Type', sep ='_')] = annotationList[[XX]][1]
            tmpEstimates[paste(Modality2_Direction,'_Type', sep ='_')] = annotationList[[XX]][2]
            tmpEstimates$InteractionType = modelAnnotations[[XX]]
            tmpEstimates <-  tmpEstimates[,sort(colnames( tmpEstimates))]
        }
                 
        return(tmpEstimates)
             
    })
             
    ## If for cytoscapeOutput, then join all the estimates together into an Edge table with annotations,
    ## And pull out row annotation data from each association model, and filter it by the edges present. 
    ## Then return a list of two tables - one for edges and one for nodes.
    if(cytoscapeOutput){
        
        edgeTable <- do.call('rbind', keptEstimates)
        rowDataList <- lapply(seq_along(keptEstimates), function(XX){
            tmpEstimates <- keptEstimates[[XX]]
            allRowData = as.data.frame(SummarizedExperiment::rowData(AssociationModelList[[XX]]))
            
            associationNames = c(paste(tmpEstimates[,5], tmpEstimates[,7], sep = "_"),
                                paste(tmpEstimates[,7], tmpEstimates[,5], sep = "_"))
            ##Filter down to associations that passed threshold.
            allRowData = allRowData[rownames(allRowData) %in% associationNames,]
            if(!'Obj2' %in% colnames(allRowData)){
            
                allRowData$Obj2 = gsub(".*_", "", rownames(allRowData))
                
            }
            #Identify row data from Obj1 and Obj2. Obj1 row data goes from 1 till the column that contains the name 'Obj1'
            #Obj2 row data goes from the column after the column named Obj1 till the column that contains the name 'Obj2'
            startCol = which(colnames(allRowData) == 'Obj1')
            endCol = which(colnames(allRowData) == 'Obj2')
            
            obj1Row = allRowData[,1:startCol, drop = FALSE]
            colnames(obj1Row) = gsub("^Obj1$", 'FeatureName', colnames(obj1Row))
            obj1Row$Modality = annotationList[[XX]][1] 
            
            obj2Row = allRowData[, c(startCol + 1):endCol, drop = FALSE]
            colnames(obj2Row) = gsub("^Obj2$", 'FeatureName', colnames(obj2Row))
            obj2Row$Modality = annotationList[[XX]][2] 
            
            plyr::rbind.fill(obj1Row, obj2Row)
        })
        
        nodeTable = dplyr::distinct(do.call(plyr::rbind.fill, rowDataList))
        
        #Re-arrange lists
        if(annotation_threshold == length(AssociationModelList)-1){
            
            type1 = grep('source__Type|Modality1__Type', colnames(keptEstimates[[1]]), value = TRUE)
            type2 = grep('target__Type|Modality2__Type', colnames(keptEstimates[[1]]), value = TRUE)
            
            cleanedDFs <- lapply(keptEstimates, function(XX){
                        
                        mod1 <- unique(XX[[type1]])
                        mod2 <- unique(XX[[type2]])
                        clean1 <- dplyr::rename(XX, !!mod1 := !!gsub("__Type","", type1), !!mod2 := !!gsub("__Type","", type2))
                        dplyr::select(clean1, !!mod1, !!mod2)
                })
            
            mergedDF2 <- cleanedDFs[[1]]
                
            for(i in 1:c(length(cleanedDFs)-1)){
                commonModality = colnames(mergedDF2)[colnames(mergedDF2) %in% colnames(cleanedDFs[[i+1]])]
                mergedDF2 <- dplyr::inner_join(mergedDF2, cleanedDFs[[i+1]], by = commonModality)
        
            }
            ## Now we can filter down the Node and Edges table to only those that appear in mergedDF2. 
            ## This will make sure that we are removing edges and nodes that aren't fully connect between associations. 
                
            #We will take the mergedDF2, which represents all the complete edges. 
            allEdges <- do.call('rbind', lapply(1:c(dim(mergedDF2)[2]-1), function(ZZ){
                
                            do.call('rbind',lapply(c(ZZ+1):dim(mergedDF2)[2], function(YY){
                                        tmp <- mergedDF2[,c(ZZ, YY)]
                                        colnames(tmp) = gsub("__Type", "", c(type1, type2))
                                        tmp
                                
                                }))
                
                }))
                
            edgeTable1 <- dplyr::semi_join(edgeTable, allEdges,  by = gsub("__Type", "", c(type1, type2)))
            colnames(allEdges) <- rev(colnames(allEdges))
            edgeTable2 <- dplyr::semi_join(edgeTable, allEdges,  by = gsub("__Type", "", c(type1, type2)))   
            
            edgeTable <- rbind(edgeTable1, edgeTable2)
                
            #Filter for only features in the edgeTable
            nodeTable <- dplyr::filter(nodeTable, FeatureName %in% unlist(allEdges))
            
                        
            #Now filter to make sure each Feature is matched to the right modality. For example, a TF and gene names can overlap. 
            fullNodeList <- rowSums(do.call('cbind', lapply(colnames(mergedDF2), function(XX){
                    
                   
                   nodeTable$FeatureName %in% mergedDF2[,XX] & nodeTable$Modality == XX
                
                }))) > 0

            nodeTable <- dplyr::distinct(nodeTable[fullNodeList,] )
            
        }
        
        if(any(duplicated(nodeTable$FeatureName))){
            
                warning("Some nodes in the node table are duplicated. Please remove duplicates manually. This occurs when different association models have contain different metadata (rowData) for the same feature. getRegulation does not know which metadata is valid and so keeps both.")
        }
            
        message(stringr::str_interp('${length(unique(nodeTable$FeatureName))} features form a network with each other.'))
        
        return(list('Edges' = edgeTable, 'Nodes' = nodeTable))
        
    } 
        
    ##If you are not returning the cytoscape output, then  merge the lists together to create one data.frame of associations.
    mergedDF <- keptEstimates[[1]]
    
    for(i in 1:c(length(keptEstimates)-1)){

        commonModality = colnames(mergedDF)[colnames(mergedDF) %in% colnames(keptEstimates[[i+1]])]
        mergedDF <- dplyr::inner_join(mergedDF, keptEstimates[[i+1]], by = commonModality)
        
    }
    totalFeatures <- lengths(lapply(unique(unlist(annotationList)), function(XX){
        
            unique(mergedDF[[XX]])
        
        }))
    message(stringr::str_interp('${sum(totalFeatures)} features form a network with each other.'))
        
    return(mergedDF)

}
                            
                            
#' @title \code{exportGeneTileLinks}
#'
#' @description \code{exportGeneTileLinksn} takes GeneTile Model and exports the significant connections to a .bedpe format or returns a data.frame.  
#' @param GeneTileAssociation A ChAI GeneTileAssociation model
#' @param filename The name of the file to be exported (should end in .bedpe if you want to view in IGV), if returnDF = FALSE
#' @param returnDF Boolean. Determines whether to write to file or return a data.frame
#' @param FDR_threshold A number between 0 and 1, used as a significance threshold for associations between genes and tiles.
#' @param filterFactors A optional list of strings that describe various factors we want to remove as sources of noise. These factors should ideally represent a strong technical artifact of somekind (like sequencing depth or other). I
#' @param backgroundThreshold A number between 0 and 1 that is used as the FDR threshold to determine whether a controlFactor is significant, and thus a given pairwise association should be removed. 
                            
#' @return a dataframe or nothing (writes to file)
#'
#'
#' @export
       
exportGeneTileLinks <- function(GeneTileAssociation, fileName,
                                returnDF = FALSE,
                                FDR_threshold = 0.1,
                                filterFactors = NULL, 
                                backgroundThreshold =0.1){
    
        tryCatch({             
            if(GeneTileAssociation@metadata$Type !=  'GeneTile'){
                
                stop('GeneTileAssociation is not actually a GeneTile Association from ChAI.')
                
            }
           },error = function(cond){
            
            stop('GeneTileAssociation is not a SummarizedExperiment from ChAI')
            
          })

    
        geneTile1 <- getEstimates(GeneTileAssociation, factor = 'exp2', FDR_threshold= FDR_threshold,
                                 filterFactors = filterFactors, backgroundThreshold = backgroundThreshold)
        geneTile2 <- getEstimates(GeneTileAssociation, factor = 'ZI_exp2', FDR_threshold= FDR_threshold,
                                 filterFactors = filterFactors, backgroundThreshold = backgroundThreshold)
        if(any(dim(geneTile1) == 0) & any(dim(geneTile2) == 0)){
            warning('No significant Gene-Tile associations found')
            return(NULL)
        }else if(any(dim(geneTile1) == 0)){
            geneTile2$Type = 'ZeroInflated'
            allGeneTiles = geneTile2
        }else if(any(dim(geneTile2) == 0)){
            geneTile1$Type = 'Continuous'
            allGeneTiles = geneTile1
        }else{
            geneTile2$Type = 'ZeroInflated'
            geneTile1$Type = 'Continuous'
            allGeneTiles = rbind(geneTile1, geneTile2)
        }

        
    
        
    
        allRowData <- SummarizedExperiment::rowData(GeneTileAssociation)
    
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
    
        write.table(allLinks, fileName, 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep ='\t')
    
        if(returnDF){
        
            return(allLinks)
           
        }
 
}
    
                            
                   
Obj1 <- Obj2 <- General <- TFs <- Type <- Genes <- Tiles <- FDR <- Estimate <- Feature <- FeatureName <- Type <- Number <- Class <- NULL