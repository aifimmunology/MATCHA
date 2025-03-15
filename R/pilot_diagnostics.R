#' @title pilotPredictions
#'
#' @description \code{pilotPredictions} extract model predictions, metadata, and adjusted measurements from the output of pilot_scRNA, pilot_scATAC, or pilot_General
#' @param modelList a list of models, output from pilotZIGLMM or pilotLMEM
#' @param  specVariable Default is NULL. If used, provide a String that matches one of the variables in your model. The experimental measurement (gene, tile, or other) will be adjusted for all other variables except for this one. 
#' @param returnList Boolean. Default is FALSE. Determines if you want a list of data.frames for each feature that was modeled, or one combined data.frame of values across features/models
#' @return list of dataframes containing metadata, and model predictions (if successful) or one combined data frame across features. 
#'
#'
#'
#' @examples
#' \dontrun{
#'   predictedDF <- plotZIModels(modelList,returnMatrix = TRUE)
#' }
#'
#' @export
#' @keywords evaluating_pilot

pilotPredictions <- function(modelList, specVariable = NULL, returnList = FALSE){

    modelCoeffs <- getPilotCoefficients(modelList)
    successfulModels = unlist(lapply(modelCoeffs, function(ZZ) { !all(is.null(ZZ))}))
    
    if(sum(successfulModels) == 0){
    
        stop('No models successfully converged')
        
    }
    
    allVariables <- rownames(modelCoeffs[[which(successfulModels)[1]]]$cond)
    allVariables <- allVariables[!grepl('ZI_|Intercept',allVariables)]
    
    metaData1 = modelList[[which(successfulModels)[1]]]$frame
    numericVariables <- allVariables[allVariables %in% colnames(metaData1)]
    remainingVariables <- allVariables[!allVariables %in% numericVariables]    
    
    if(!is.null(specVariable)){
        if(!specVariable %in% allVariables){

            stop('specVariable not found in model. Make sure spelling is correct, and double check categorical variable names')

        }                        
        adjustVariables <- allVariables[! allVariables %in% specVariable]   
    }else{
        
        adjustVariable = allVariables
        
    }

    
    if(any(grepl('list|List',class(modelList)[1]))){
        
      pList <- lapply(seq_along(modelList), function(i){
          
          if(methods::is(modelList[[i]], 'glmmTMB')){
            #Extract the original values
            df <- as.data.frame(modelList[[i]]$frame)
            feature <- names(modelList)[i]
              
            ## Log transform the FragmentNumbers so as to stabilize the model. But only if FragmentCounts is in the model. Same for CellCounts.
            if(any(colnames(df ) %in% c('FragmentCounts'))){
                df$rawFragmentCounts = df$FragmentCounts
                df$FragmentCounts <- log10(df$FragmentCounts)
            }
            if(any(colnames(df ) %in% c('CellCounts'))){
                df$rawCellCounts = df$CellCounts
                df$CellCounts <- log10(df $CellCounts)
            }

            #Create a metadata column for each categorical variable, one-hot encoding them. 
            if(length(remainingVariables) > 1){

                nextVariables <- lapply(remainingVariables, function(x) matchCategorical(df, x))
                newMetaData = do.call('cbind', nextVariables)
                colnames(newMetaData) <- remainingVariables

            }else if(length(remainingVariables) > 0){
                newMetaData <- data.frame(newVar = matchCategorical(df, remainingVariables))
                colnames(newMetaData) <- remainingVariables
            }else{
                newMetaData <- NULL
            }


            df = cbind(df, newMetaData)  
            df$Measurement = feature
                      
            ## Extract coefficients, or simply return DF with blank predictions if failed to converge.
            tryCatch({
              #Extract model coefficients
              sum_fit = summary(modelList[[i]])
              coefs = sum_fit$coefficients$cond[,1, drop= FALSE]
              
              allAdjusts = do.call('cbind',lapply(allVariables[! allVariables %in% specVariable], function(x){

                            as.data.frame(coefs[x,]*as.numeric(df[,x]))

                            }))
                
                if(!is.null(specVariable)){

                    predictVariables <- as.matrix(df[, specVariable, drop = FALSE]) %*% as.matrix(coefs[specVariable,, drop = FALSE])
                    exp_adjusted = df$exp - rowSums(allAdjusts)
                    Prediction = coefs[grepl('Intercept',rownames(coefs)),] + predictVariables
                    df <- cbind(data.frame('Prediction' = Prediction, 'exp_adjusted' = exp_adjusted, 
                                          'Estimate' = rep(unlist(coefs[specVariable,1]), dim(df)[1]),
                                          'PValue' = rep(unlist(coefs[specVariable,1]), dim(df)[1])), df)


                }else{

                    Prediction = coefs[grepl('Intercept',rownames(coefs)),] + rowSums(allAdjusts)
                    df <- cbind(data.frame('Prediction' = Prediction), df)

                 }
            
                return(df)
           
                },error=function(e){
                
                     if(!is.null(specVariable)){
                     
                         df <- cbind(data.frame('Prediction' = rep(NA, dim(df)[1]), 'exp_adjusted' = rep(NA, dim(df)[1]), 
                                          'Estimate' = rep(NA, dim(df)[1]), 'PValue' = rep(NA, dim(df)[1])), df)
                     
                    }else{

                         df <- cbind(data.frame('Prediction' = rep(NA, dim(df)[1])), df)
                    }
                    
                })
                                        
            }else{
                #If the model completely failed, the return the data frame. 
                df <- as.data.frame(modelList[[i]]$Data)
                feature <- modelList[[i]]$Measurement
                if(!is.null(specVariable)){
                     
                     df <- cbind(data.frame('Prediction' = rep(NA, dim(df)[1]), 'exp_adjusted' = rep(NA, dim(df)[1]),  
                                            'Estimate' = rep(NA, dim(df)[1]), 'PValue' = rep(NA, dim(df)[1])), df)
                     
                }else{
                     
                     df <- cbind(data.frame('Prediction' = rep(NA, dim(df)[1]),  
                                            'Estimate' = rep(NA, dim(df)[1]), 'PValue' = rep(NA, dim(df)[1])), df)
                }
                
            }
                                     
            #return the data.frame of data, with predictions if the model converged. 
            return(df)
                                        
        })
                
      #return the list of data.frames
    if(returnList){
    
        return(pList)
    
    }else{
         
        return(do.call('rbind', pList))
        
    }
                
    }else{
      stop('modelList type not recognized.')
    }
  }



#' @title getPilotCoefficients
#'
#' @description \code{getPilotCoefficients} Attempts to pull coefficients from a list of models from pilotZIGLMM or pilotLMEM. 
#'    Returns a list of either the coefficients for each
#'    model, or the error generated when attempting to get coefficients. 
#' @param pilotModelList A list of models, output from any pilot modeling functions. 
#' @param returnSummary A boolean, to determine whether to return the individual summaries of each model, or
#'         to compile those results into a single data.frame for browsing. Default is TRUE (returns individual summaries) 
#' @param paramInterest A string, only used when returnSummary = FALSE. The string describes which 
#'       coefficient parameter you want summarized into a data.frame. 
#'       Options include 'estimate', 'std_error', 'z_value', and 'p_value'
#' @return modelList a list of outputs from glmmTMB
#'
#'
#' @export
#' @keywords evaluating_pilot

getPilotCoefficients <- function(pilotModelList, returnSummary = TRUE, paramInterest = 'p_value'){

   coeffList <-  lapply(pilotModelList, function(x) {
                    if(length(x) == 7){
                        tryCatch({
                            summary(x)$coefficients
                        },
                        error=function(e){
                            list(e, x$frame)
                        })
                    }else{
                        x[c(1,3)]
                    }
                  
                  })
   if(returnSummary){
       return(coeffList)
    }
    
    options1 = c(1:4)
    names(options1) = c('estimate', 'std_error', 'z_value', 'p_value')
       
    if(!tolower(paramInterest) %in% names(options1)){
        
        stop('paramInterest not recognized. Must be estimate, std_error, z_value, or p_value')
        
    }
    #Identify which column index to pull out
    specIndex = which(names(options1) %in% tolower(paramInterest))
        #interate over all models. Pull out the parameter of interest and merge into a data.frame
    coeffList2 <- lapply(coeffList, function(XX){
       do.call('c', lapply(c(1:2), function(Z){
               if(!is.null(XX[[Z]]) & methods::is(XX[[Z]], 'matrix')){
                   df <- as.data.frame(XX[[Z]])[,specIndex]
                   if(Z == 2){
                       names(df) = paste('ZI', rownames(XX[[Z]]), sep = "_")
                   }else if(Z==1){
                        names(df) = rownames(XX[[Z]])
                   }
                   return(df)
                }else{
                   return(NULL)
                }
           }))
    })
    coeffDF <- suppressWarnings(do.call('rbind', coeffList2))
    
   return(as.data.frame(coeffDF))

}


#' @title getConvergenceRate
#'
#' @description \code{getConvergenceRate} Look at the model output list for any modeling function (whether unimodal or associative),
#'                      And get the rate at which the glmmTMB successfully converges on coefficients. 
#' @param modelList a list of models, output from pilotZIGLMM or pilotLMEM
#' @return nothing. Prints the convergence rate out as a message. 
#'
#'
#'
#' @examples
#' \dontrun{
#'   getConvergenceRate(modelList)
#' }
#'
#' @export
#' @keywords evaluating_pilot

getConvergenceRate <- function(modelList){
    
    if(all(unlist(lapply(modelList, function(XX)methods::is(XX, 'list')))) & 
        !('error' %in% names(modelList[[1]]))){
        
       tmp1 = lapply(names(modelList), function(XX){
               print(XX)
               getConvergenceRate(modelList[[XX]])
           })
           
    }else{ 
        modelOuts <- getPilotCoefficients(modelList)
        ## Identify types of model outputs. If it's a list of length 3, then it'll be 
        modelOutLengths <- lengths(modelOuts)

        mainFailure = sum(modelOutLengths !=3)
        otherFailure <- sum(unlist(lapply(which(modelOutLengths == 3), function(XX){
                if(dim(modelOuts[[XX]]$cond)[1] > 1){
                    modelOutTmp <- modelOuts[[XX]]$cond[,2]
                    all(is.na(modelOutTmp[names(modelOutTmp) != '(Intercept)']))
                }else{
                    is.na(modelOuts[[XX]]$cond[,2])
                }
            })))

        totalNumber = length(modelList)

        successNumber = totalNumber-sum(mainFailure + otherFailure)

        message('Convergence Rate: ', successNumber/totalNumber*100, '%')
    }

}


#' @title plotPilotFits
#'
#' @description \code{plotPilotFits} Runs DHARMa::simulateResiduals for zero-inflated models. 
#' @param modelList A list of models from any ChAI pilot modeling function (whether on a single modality or the association between two modalities)
#' @param fileName Name of the pdf to export. The function will add .pdf to the end of this string.
#' @return NULL Save a pdf under the given name with residuals plotted. 
#'
#'
#'
#' @examples
#' \dontrun{
#'   pilot_performance(modelObj)
#' }
#'
#' @export
#' @keywords evaluating_pilot
#'

plotPilotFits <- function(modelList, fileName){

    pdfName2 <- paste(fileName, '.pdf', sep ='')

    grDevices::pdf(pdfName2)
    
    if(all(unlist(lapply(modelList, function(XX)methods::is(XX, 'list')))) & 
        !('error' %in% names(modelList[[1]]))){
        
       lapply(modelList, plotPilotFits)
           
    }else{

        for(i in 1:length(modelList)){

             tryCatch({
                    print(DHARMa::simulateResiduals(modelList[[i]], plot= T))
                }, error = function(e){
                    return(e)
                })
        }
    }

   
    grDevices::dev.off()
    
    
    return(NULL)

}





#' @title testPilotFits 
#'
#' @description \code{testPilotFits } Runs DHARMa::simulateResiduals for zero-inflated models, and runs a variety of tests on each model to evaluate model fit.
#' @param modelList A list of models from any ChAI pilot modeling function (whether on a single modality or the association between two modalities)
#' @param numCores Number of cores to parallelize over.
#' @return NULL. Save a pdf under the given name with residuals plotted. 
#'
#'
#'
#' @examples
#' \dontrun{
#'   pilot_performance(modelObj)
#' }
#'
#' @export
#' @keywords evaluating_pilot

testPilotFits <- function(modelList, numCores = 1){

    if(all(unlist(lapply(modelList, function(XX)methods::is(XX, 'list')))) & 
        !('error' %in% names(modelList[[1]]))){
        
       allRes <- lapply(modelList, testPilotFits, numCores = numCores)
       names(allRes) = names(modelList)
       return(allRes)
           
   }
    
    if(numCores <= 1){
        cl = NULL
    }else{

        cl <- parallel::makeCluster(numCores)
        parallel::clusterEvalQ(cl, {
            library('DHARMa')
          })
        
    }
  # Make your clusters for efficient parallelization
   testFits <- function(i){
    
        tryCatch({
            testFit <- DHARMa::simulateResiduals(i, plot= F)
            df <- data.frame(
                KSTest_PValue = stats::ks.test(testFit$observedResponse, testFit$simulatedResponse)$p.value,
                Dispersion_PValue = DHARMa::testDispersion(testFit, plot = F)$p.value,
                Outlier_PValue = DHARMa::testOutliers(testFit, plot = F)$p.value,
                Uniformity_PValue = DHARMa::testUniformity(testFit, plot = F)$p.value,
                Quantile_PValue = DHARMa::testQuantiles(testFit, plot = F)$p.value,
                ZeroInflation_PValue = DHARMa::testZeroInflation(testFit, plot = F)$p.value
            )
        }, error = function(e){
           NULL
        })
        if(methods::is(df, 'data.frame')){
            return(df)
        }else{
            
            return(NULL)
            
        }
        
    }
  
    residualFits <- pbapply::pblapply(cl = cl, X= modelList, testFits)
    residualFits <- do.call('rbind', residualFits)
    residualFits$Feature = rownames(residualFits)
    
    
    if(numCores > 1){
        parallel::stopCluster(cl)
    }
    
    return(residualFits)

}

#' @title drop1Analysis
#'
#' @description \code{drop1Analysis} Runs drop1 analysis across features in a modelList
#' @param modelList A list of models from any ChAI pilot modeling function (whether on a single modality or the association between two modalities). It should not be zero-inflated.
#' @param numCores The number of cores to paralelize over
#' @return A dataframe of drop1Analysis results over all features, where rows are features. Failed models are not included
#'
#'
#'
#' @examples
#' \dontrun{
#'   drop1Analysis(modelObj)
#' }
#'
#' @export
#' @keywords evaluating_pilot

drop1Analysis <- function(modelList, numCores = 1){
    
    message('Running drop1Analysis, assuming that the models only involve a continuous component.')
    
    if(all(unlist(lapply(modelList, function(XX)methods::is(XX, 'list')))) & 
        !('error' %in% names(modelList[[1]]))){
        
       allRes <- lapply(modelList, drop1Analysis, numCores = numCores)
       names(allRes) = names(modelList)
       return(allRes)
           
   }
    
    if(numCores <= 1){
        cl = NULL
    }else{

        cl <- parallel::makeCluster(numCores)
       
    }
        browser()
    model_drops <-  pbapply::pblapply(cl = cl, X= modelList, testDrop)
    
    modelDF <- do.call('rbind', model_drops)
    modelDF$Feature = names(modelList)[!unlist(lapply(model_drops, function(XX) {all(is.null(XX))}))]
    
    if(numCores > 1){
        parallel::stopCluster(cl)
    }
    
    return(modelDF)

}

#' @title internal drop1Analysis function
#'
#' @description \code{testDrop} Runs drop1 analysis one one model
#' @param res A glmmTMB object
#' @return A dataframe of drop1Analysis results for one model (i.e. feature)
#'
#'
#'
#' @examples
#' \dontrun{
#'   drop1Analysis(modelObj)
#' }
#'
#' @noRd
testDrop <- function(res){
         df <- tryCatch({
            drop1Test <- as.data.frame(stats::drop1(res, test = 'Chisq'))
            drop1Test$Factor = rownames(drop1Test)
            drop1Test$Factor[1] = 'None'
            colnames(drop1Test) <- c('Df', 'AIC','LRT', 'P_Value', 'Factor')
            tidyr::pivot_wider(drop1Test[,c('AIC', 'P_Value','Factor')], 
                            names_from = 'Factor', 
                            values_from = c('AIC', 'P_Value'))
        }, error = function(e){
           NULL
           
        })
        return(df)
    }
          
#' @title summarizeResFits
#'
#' @description \code{summarizeResFits} Summarized the output of testPilotFits, allowing a quick summary of how many models (as a percentage) failed which tests. 
#' @param summaryOfFits A data frame, generated by testPilotFits
#' @param threshold A p-value threshold used for the various tests, to check whether the model for a feature failed a test. The standard is p = 0.05.
#' @param uncorrected A boolean, indicating whether to use the nominal p-value for each statistical tests, or to convert those p-values to a corrected p-value before calculating the percentage of models that failed a given test. Default is TRUE. If FALSE, then corrected p-values are used. 
#' @param correctionType A string, documenting which type of p-value adjustment to use, if uncorrected = FALSE. See R documentation on p.adjust for options. 
#' @param na.rm Some tests may generate NAs. This removes those NAs when calculating the percent failed. Default is FALSE. 
#' @param verboseOutput A boolean flag, that determines whether to print the results as a percentage, or return a named vector with numeric values. 
#' @return A named lists, showing what percentage of models failed a given test. 
#'
#'
#'
#' @examples
#' \dontrun{
#'   summarizeResFits(model_fit_data_frame)
#' }
#'
#' @export
#' @keywords evaluating_pilot

summarizeResFits <- function(summaryOfFits, 
                             threshold = 0.05,
                             uncorrected = TRUE,
                             correctionType = 'fdr',
                             verboseOutput = TRUE,
                             na.rm = FALSE){
    
    if(!methods::is(summaryOfFits, 'data.frame')){
        
       allRes <- do.call('cbind',lapply(summaryOfFits, function(ZZ){
               summarizeResFits(ZZ, 
                             threshold = threshold,
                             uncorrected = uncorrected,
                             correctionType = correctionType,
                             verboseOutput = verboseOutput,
                             na.rm = na.rm)
           }))
    
       colnames(allRes) = paste(names(summaryOfFits),'Failed',sep = "_")
       return(allRes)
           
   }
    
    
    if(!methods::is(summaryOfFits, 'data.frame')){
        stop('summaryOfFits is not a data.frame')
    }
    reqCols = c("KSTest_PValue", "Dispersion_PValue", "Outlier_PValue", 
                "Uniformity_PValue", "Quantile_PValue","ZeroInflation_PValue", "Feature")
    if(!all(reqCols == colnames(summaryOfFits))){
        stop('summaryOfFits does match the format (column names) that comes from testPilotFits().')
    }
       
    if(uncorrected){
        res1 <- unlist(apply(summaryOfFits[,c(1:6)], 2, function(x) 
                        sum(x <threshold, na.rm = na.rm)/length(x)))
    }else {
        
        res1 <- unlist(apply(summaryOfFits[,c(1:6)], 2, function(x){
                        x2 <- stats::p.adjust(x, method = correctionType)
                        
                        sum(x2 <threshold, na.rm = na.rm)/length(x2)
            }))
        
    }
                             
    if(verboseOutput){
        
        res2 <- paste0(formatC(100 * res1, format = "f", digits = 1), "%")
        names(res2) <- gsub("_PValue", "", names(res1))
        res2 <- as.data.frame(res2)
        colnames(res2) <- 'PercentFailed'
        return(res2)
        
    }else{
        names(res1) <- gsub("_PValue", "", names(res1))
        return(res1)
    }

}
