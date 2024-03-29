# ----------------------------------------------------------------------- #
# Geographic Patterns in U.S. Lung Cancer Mortality and Cigarette Smoking
# ----------------------------------------------------------------------- #
#
# Created by: Ian Buller, Ph.D., M.A. (GitHub: @idblr)
# Created on: January 15, 2022
#
# Most recently modified by: @idblr
# Most recently modified on: May 17, 2022
#
# Notes:
# A) Code to prepare the data for the analysis
# B) Uses the soon-to-depreciate 'sp' package (because of 'spdep' package dependency for Lee's L calculation)
# C) 05/17/2022: Include 'spatialreg' package for symmetric matrices; remove "raster", "rgeos", and "rgdal" packages
# D) 05/17/2022: Consistent variable name based on merged data object
# E) 05/17/2022: Use 'sf' objects instead of 'sp' objects when available
# F) 05/17/2022: Cropping by coastline data no longer necessary with "cb = TRUE" for 'tigris' data
# G) 05/17/2022: Rename NAs as "Suppressed" for LeeFDR results
# ----------------------------------------------------------------------- #

cat("\nPreparing data for Local Lee's L statistic\n")

############
# PACKAGES #
############

loadedPackages <- c("boot", "cowplot", "dplyr", "ggplot2", "grDevices", "gtools",
                    "Matrix", "sf", "sp", "spatialreg","spdep", "stats", "stringr",
                    "tigris", "utils")
suppressMessages(invisible(lapply(loadedPackages, require, character.only = TRUE)))
options(tigris_use_cache = T)

############
# SETTINGS #
############

# RNG seed 
set.seed(15335) # reproducibility with manuscript
# # uncomment for new random seed
# initial_seed <- as.integer(Sys.time())
# the_seed <- initial_seed %% 100000 # take the trailing five digits of the initial seed
# set.seed(the_seed)

# Alpha level
alpha <- 0.05  # for a 95% confidence interval

####################
# CUSTOM FUNCTIONS #
####################

source("code/functions.R")

####################
# DATA IMPORTATION #
####################

# Path for county-level lung cancer mortality rate and smoking prevalences
dat_path <- "data/lung_cancer_mortality_and_smoking_prevalence.csv"

# Generate one CSV files
## 'data.frame': 3113 obs. of  8 variables:
# $ FIPS                : chr  [FIPS codes for each county]
# $ Female_Rate         : chr  [Female Lung Cancer Mortality Rate]
# $ Male_Rate           : chr  [Male Lung Cancer Mortality Rate]
# $ ESF_Pct_1997_2003   : num  [Female Ever Smoking Prevalence]
# $ ESM_Pct_1997_2003   : chr  [Male Ever Smoking Prevalence]
# $ CSF_Pct_1997_2003   : chr  [Female Current Smoking Prevalence]
# $ CSM_Pct_1997_2003   : num  [Male Current Smoking Prevalence]

# import the .csv file "Lung_Cancer_Mortality_and_Smoking Prevalence.csv" 
alcove_dat <- utils::read.csv(dat_path)
#utils::str(alcove_dat)

## Summary
# n = 3,113 counties
# 75 variables

# Fix FIPS in dataset to include leading 0 for states with STATE FIPS < 10
alcove_dat$FIPS <- stringr::str_pad(alcove_dat$FIPS, 5, "left", "0")

##################
# SPATIAL WINDOW #
##################

# US State Shapefile
not_l48 <- c("Commonwealth of the Northern Mariana Islands", "Guam", "American Samoa",
             "Hawaii", "Alaska", "Puerto Rico", "United States Virgin Islands")
shp_state <- suppressMessages(tigris::states(year = 2018, class = "sf", cb = TRUE))
shp_state_name <- shp_state$NAME[shp_state$NAME %notin% not_l48]

# Lower 48 US State Shapefile (as 'sf' object for plotting)
shp_l48 <- shp_state[shp_state$NAME %notin% not_l48, ]
proj_l48 <-  sf::st_transform(shp_l48, crs = 4326)

# Lower 48 US County Shapefile (as 'sp' object for Lee's L statistic)
shp_US <- suppressMessages(tigris::counties(state = shp_state_name, year = 2018, class = "sp", cb = TRUE))
proj_US <- sp::spTransform(shp_US, sp::CRS("EPSG:4326"))

###################
# SPATIAL LINKAGE #
###################

# Recode Oglala Lakota County, South Dakota as Shannon County, South Dakota
proj_US[proj_US$NAME == "Oglala Lakota", "GEOID"] <- "46113"

# Merge rate data with county shapefile
alcove_geo <- sp::merge(proj_US, alcove_dat, 
                        by.x = "GEOID", by.y = "FIPS",
                        all.x = TRUE, all.y = TRUE)

######################
# SPATIAL PROJECTION #
######################

kansas_epsg <- "EPSG:26978" # US geographical center near Kansas
alcove_proj <- sp::spTransform(alcove_geo, sp::CRS(kansas_epsg))

############################
# LEE's (2001) L STATISTIC #
############################
# a bivariate spatial association measure called the L statistic that integrates Pearson's r and Moran's I statistics
# account for both spatial and in-situ correlations
# an improvement upon bivariate LISA
# https://doi.org/10.1007/s101090100064

cat("\nCalculating Local Lee's L for Female Lung Cancer Mortality Rate by Female Current Smoking Prevalence... 1 of 4")
LeeCSF <- LeeL(dat = alcove_proj, x = "CSF_Pct_1997_2003", y = "Female_Rate", label = "LCSF", numsim = 100000)
cat("\nCalculating Local Lee's L for Female Lung Cancer Mortality Rate by Female Ever Smoking Prevalence... 2 of 4")
LeeESF <- LeeL(dat = alcove_proj, x = "ESF_Pct_1997_2003", y = "Female_Rate", label = "LESF", numsim = 100000)
cat("\nCalculating Local Lee's L for Male Lung Cancer Mortality Rate by Male Current Smoking Prevalence... 3 of 4")
LeeCSM <- LeeL(dat = alcove_proj, x = "CSM_Pct_1997_2003", y = "Male_Rate", label = "LCSM", numsim = 100000)
cat("\nCalculating Local Lee's L for Male Lung Cancer Mortality Rate by Male Ever Smoking Prevalence... 4 of 4")
LeeESM <- LeeL(dat = alcove_proj, x = "ESM_Pct_1997_2003", y = "Male_Rate", label = "LESM", numsim = 100000)
cat("\nCleaning-up")

Lee <- merge(LeeCSF$out, LeeESF$out, by = "GEOID")
Lee <- merge(Lee, LeeCSM$out, by = "GEOID")
Lee <- merge(Lee, LeeESM$out, by = "GEOID")
Lee <- suppressMessages(sf::st_as_sf(Lee)) # Convert to an 'sf' object

FDRpvals <- data.frame(x = c("LCSF", "LESF", "LCSM", "LESM"),
                       y = c(LeeCSF$FDRpval, LeeESF$FDRpval, LeeCSM$FDRpval, LeeESM$FDRpval))
names(FDRpvals) <- c("Comparison", "FDR_pvals")

# Rename NAs as "Suppressed
Lee$LESM_LeeFDR[is.na(Lee$LESM_LeeFDR)] <-  "Suppressed"
Lee$LESF_LeeFDR[is.na(Lee$LESF_LeeFDR)] <-  "Suppressed"
Lee$LCSM_LeeFDR[is.na(Lee$LCSM_LeeFDR)] <-  "Suppressed"
Lee$LCSF_LeeFDR[is.na(Lee$LCSF_LeeFDR)] <-  "Suppressed"
  
###########
# CLEANUP #
###########

rm(list = setdiff(ls(), c("Lee", "FDRpvals", "proj_l48")))

# ----------------------------- END OF CODE ----------------------------- #
