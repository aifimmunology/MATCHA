# MATCHA: Modeling Associations for all The Creative, Heterogeneous Assays

MATCHA explicitly moves beyond statistical modeling of any one type of biological data types (scRNAseq, scATACseq, Olink, flow etc...). Instead, it provides both generic and specialized functions for :

* High-throughput modeling on high dimensional data types: 
* High-throughput modeling of relationships between features of different data types (but from the same biological sample)
  * Examples:
    * Protein A (Serum Proteome) vs. Cell Type B Frequency (Flow Cytometr)
    * Region C Accessibility (scATAC-seq) vs. Gene D Expression (scRNA-seq)
* Visualizing, evaluating, exporting, and integrating model results

Furthermore, all modeling is based on glmmTMB, allowing for:

* Cross-sectional analysis
* Longitudinal analysis
* Zero-inflated modeling
* A variety of data distributions (see glmmTMB families)

As a result, the user can test for complex relationships, and generate statistically-informed causal hypotheses for their biological questions (gene regulatory networks, signaling cascades, drivers of proliferation, or disease progression, etc... ). Almost any data type can be imported to matcha, which has been vetted with:

* scRNA-seq (Seurat)
* scATAC-seq (MOCHA, directly compatible)
* ChromVAR (via MOCHA)
* proteomics (olink)
* Cell frequencies/abundances

See overview for more details, as well as the [Tips & Tricks article]()

------------------------------------------------------------------------

### Table of Contents

-   [Installation](#installation)
-   [Overview](#overview)
-   [Contact](#contact)


-----------------------------------------------------------------------

## <a name="installation"></a> Installation

Devtools installation requires CMake, which can be installed on Ubuntu-based systems via:

    sudo apt-get -y install cmake

For other systems, look here: https://cgold.readthedocs.io/en/latest/first-step/installation.html

Install from GitHub:

    devtools::install_github("aifimmunology/matcha")

Or install a specific development branch from GitHub:

    devtools::install_github("aifimmunology/matcha", ref = "your_branch_name")

The development branch has the latest stable version of matcha.

    devtools::install_github("aifimmunology/matcha", ref = "development")
    
## <a name="overview"></a> Overview

The package includes a set of funcitons for:

* Data import (scRNA, ChromVAR, or any other modality. MOCHA's scATAC-seq is directly compatible)
    * makePseudobulkRNA, makeChromVAR, importGeneralModality
* Testing out formulas ("pilotting" a given formula to determine which variables should be modeled)
    * pilot_scRNA/scATAC/ChromVAR/General, pilot_GeneTile/scRNA_Associations/scATAC_Assocations/ChromVAR_Assocations/General_Associations
* Evaluating model fits for all pilot models (built off of Dharma's functions)
    * getConverageRate, getPilotCoefficients, testPilotFits, plotPilotFits
* Modeling an individual data type
    * model_scRNA/scATAC/ChromVAR/General
* Modeling assocations between features in different data types
    * GeneTile_Associations, scATAC_Associations, scRNA_Associations, ChromVAR_Associations, General_Associations
* Getting & transforming Data Objects
    * subsetmatcha, mergematcha, flattenmatcha, transformmatcha, renameCellTypes
    * getDetectionRates, getAdditionalData, getPopulationMetaData
* Getting & using model information
    * getEstimates, modelPredictions, associatedPredictions
    * getResiduals, getModelFactors, getModelValues
* Network creation & export
    * getRegulation, exportGeneTileLinks


![Workflow Diagram](man/figures/matchaworkflow.png)
