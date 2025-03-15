### Test for pilot models

test_that("Pilot models for scATAC works.", {
    
 capture.output(
  p1 <- pilot_scATAC(ChAI:::scATAC_Object,
                    cellPopulation = 'CD16 Mono',
                    continuousFormula =  exp ~ Age + Sex + Group + (1|PTID),
                    ziFormula = ~ 0 + FragNumber + Group,
                    zi_threshold = 0,
                    verbose = FALSE,
                    numCores = 3,
                    featureSample = 20),
      type = "message"
    )
    
  expect_snapshot(
      p1
    )
})


tmp <- readRDS('../COVID_scATAC_Manuscript/ChAI_Data_scRNA.rds')
scRNA_Object <- subsetChAI(tmp, subsetBy = 'PTID_visit', groupList = lateDiffDF$PTID_visit)
usethis::use_data(scRNA_Object)

tmp2 <- readRDS('../COVID_scATAC_Manuscript/Olink_SummarizedExperiment.rds')
Olink_Object <- subsetChAI(tmp2, subsetBy = 'PTID_visit', groupList = lateDiffDF$PTID_visit)
usethis::use_data(Olink_Object)

test_that("Pilot models for scRNA works.", {
    
 capture.output(
  p1 <- pilot_scRNA(ChAI:::scRNA_Object,
                    cellPopulation = 'CD16 Mono',
                    modelFormula =  exp ~ Age + Sex + Group + (1|PTID),
                    verbose = FALSE,
                    numCores = 3,
                    featureSample = 20),
      type = "message"
    )
    
  expect_snapshot(
      p1
    )
})


test_that("Pilot models for General works.", {
    
 capture.output(
  p1 <- pilot_scRNA(ChAI:::Olink_Object,
                    cellPopulation = 'CD16 Mono',
                    modelFormula =  exp ~ Age + Sex + Group + (1|PTID),
                    verbose = FALSE,
                    numCores = 3,
                    featureSample = 20),
      type = "message"
    )
    
  expect_snapshot(
      p1
    )
})

test_that("Pilot models for ChromVAR works.", {
    
 capture.output(
  p1 <- pilot_ChromVAR(ChAI:::Olink_Object,
                    cellPopulation = 'CD16 Mono',
                    modelFormula =  exp ~ Age + Sex + Group + (1|PTID),
                    verbose = FALSE,
                    numCores = 3,
                    featureSample = 20),
      type = "message"
    )
    
  expect_snapshot(
      p1
    )
})

