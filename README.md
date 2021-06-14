Erickson, Kelley D. & Adam B. Smith Accounting for imperfect detection in data from museums and herbaria when modeling species distributions: Combining and contrasting data-level versus model-level bias correction. Ecography 

Preprint available at https://www.biorxiv.org/content/10.1101/2021.01.06.425644v1 


**Code**

- **01_Download_Data.R**: Script for downloading the raw data (GBIF occurrences, GADM shapefile, and 2010 census data). Creates: `tracheophytes_coords`, `tracheophytes_nocoords`, `florida_state`, `florida` `census_2010`  
- **02_Process_Data.Rmd**: Script for processing the raw data. Creates: `surv_cov3.RData` and `unit_cov.RData`
- **03_Assemble_Simple_Model_Folds.R**: Generates cross-validation folds for the simple models
- **04A_Crossvalidation_FullModel.R**: Runs cross-validation for the full model using the full data set or the filtered data set. 
- **04B_Crossvalidation_SimpleModel.R**: Runs cross-validation for the simple model using the full data set or the filtered data set. 
- **05_Figures.R**: Code for generating the figures in the manuscript and supplements. 


**Data files **
- **correctedDates2.csv**: A csv file of corrections for collector name or year that were identified by visually examining specimen sheets virtually. 3 columns: gbifID (the unique identifier assigned by GBIF to an occurrence record), correctedYear (the corrected year), and correctedCollector (a corrected collector name). NA values indicate that it was not possible to identify the correct year or collector from the herbarium record. 
- **flagged_coordinate_data.csv**: A csv file of occurrence records where the county name in the county field (from the GBIF occurrence data) did not match the geographically assigned county (from data processing). 4 columns: gbifID(the unique identifier assigned by GBIF to each occurrence record), county(the county name as recorded in GBIF), geo_county(the name of the GADM county in which the occurrence was determined to fall based on coordinates), FLAG("FLAG" means that the county and geo_county do not match and this occurrence record should be removed from the analysis, "" means that the county and geo_county do match and this occurrence record should not be removed from analysis
- **incorrectYearsCollectors.csv**: CSV file created by examining occurrences with collection years that were potentially outside of the realm of realistic collection years for a particular collector (i.e. after the known death date of a collector, or several years before or after other collections by a collector) 4 columns: GBIF_ID (the unique identifier assigned by GBIF to each occurrence record), CorrectedYear(either: the corrected year for a given occurrence record based on verification by looking at the herbarium sheet, "" to indicate that no correction is needed, or FLAG to indicate that an occurrence record should not be used in the analysis because the date was not verifiable), CorrectedCollector(either: the corrected collector name for a given occurrence record based on verification by looking at the herbarium sheet, "" to indicate that no correction was needed, or FLAG to indicate that an occurrence record should be removed from the analysis since the collector name could not be verified), Notes (links to herbarium specimen sheets)
- **listCollectors_refined.csv**:A lookup table that includes collector names that were corrected and cleaned using openRefine. 3 columns: Column: numeric rownames for each occurrence, gbifID(the unique identifier assigned by GBIF to each occurrence record), recordedBy (the corrected collector names)
- **MANIND\**
  + **collectors_of_species.RData**
  + **surv_cov.RData**
  + **surv_cov3.RData**
- **METTOX\**
- **RHUCOP\**
- **SCHTER\**
- **TOXPUB\**
- **TOXRAD\**
- **TOXVER\**
-**unit_cov2.RData**