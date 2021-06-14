Erickson, Kelley D. & Adam B. Smith Accounting for imperfect detection in data from museums and herbaria when modeling species distributions: Combining and contrasting data-level versus model-level bias correction. Ecography 

Preprint available at https://www.biorxiv.org/content/10.1101/2021.01.06.425644v1 


**Code**

- **01_Download_Data.R**: Script for downloading the raw data (GBIF occurrences, GADM shapefile, and 2010 census data). Creates: `tracheophytes_coords`, `tracheophytes_nocoords`, `florida_state`, `florida` `census_2010`  
- **02_Process_Data.Rmd**: Script for processing the raw data. Creates: `surv_cov3.RData` and `unit_cov.RData`
- **03_Assemble_Simple_Model_Folds.R**: Generates cross-validation folds for the simple models
- **04A_Crossvalidation_FullModel.R**: Runs cross-validation for the full model using the full data set or the filtered data set. 
- **04B_Crossvalidation_SimpleModel.R**: Runs cross-validation for the simple model using the full data set or the filtered data set. 
- **05_Figures.R**: Code for generating the figures in the manuscript and supplements. 


