# ChAI 1.0.0

ChAI now includes a large number of new functions for data handling, model evaluation, network construction, data export, browsing, and more.

* Added testPilotFits(), plotPilotFits(), drop1Analysis(), and summarizeResFits() for evaluating model performance
* ChAI now explicitly tracks data and model types to ensure compatibility. 
    * All functions now test for ChAI data object types (scRNA, scATAC, ChromVAR, or General), and check for compatibility with the functions being used.
    * All functions that operate on ChAI model objects also take into account the type of model involved, allowing for data-type specific transformations. 
* ChAI now has explicity functions for ChromVAR modeling, rather than relying on the general modeling function. This allows technical impacts, like FragNumber to be included in the modeling formula.
* A new function (getRegulation) building networks from 3 or more associative models (minimum of 3 modalities). It can export node and edge tables for network analysis in Cytoscape or other software packages as well. 
* A new function (exportGeneTileLinks) can now export Gene-Tile associations into a .bedpe file for browsing in IGV viewer or other softwares. 
* A name change was made - functions for using model outputs and data are now labeled modelPredictions (unimodal models) and associatedPredictions (association models), to avoid confusion with other functions. These two functions allow the user to plot model predictions, and correct data for sources of technical or biological influence. 
* A number of object utils have also been added that allow the user to
    * access model factors (getModelFactors)
    * rename celltypes in a scRNA, scATAC, or ChromVAR data object
    * access population-level metadata (i.e. CellCounts, average reads per cell, etc.. ) for scATAC, scRNA, and ChromVAR (take from scATAC pseudobulk)
    * access the detectionRate for genes across cell populations (scRNA) - see getDetectionRates()


# ChAI 0.0.1

* Added a `NEWS.md` file to track changes to the package.
