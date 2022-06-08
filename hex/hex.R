# ------------------------------------------------------------------------------ #
# Hexsticker for the GitHub Repository idblr/geo_US_lung_cancer_and_smoking
# ------------------------------------------------------------------------------ #
#
# Created by: Ian Buller, Ph.D., M.A. (GitHub: @idblr)
# Created on: May 17, 2022
#
# Recently modified by: @idblr
# Recently modified on: June 08, 2022
#
# Notes:
# A) Uses the "hexSticker" package
# B) Modified image from a companion manuscript figure
# D) Hexsticker for the GitHub Repository https://github.com/idblr/geo_US_lung_cancer_and_smoking
# ------------------------------------------------------------------------------ #

####################
# DATA IMPORTATION #
####################

# Use the code found in 'functions.R' and 'preparation.R' files
## Loads three objects
### A) 'Lee' an 'sf' object of the U.S. county-level lung cancer mortality rates, smoking prevalences, and Lee's L values
### B) 'FDRpvals' a 'data.frame' object of significant p-value cutoff for each bivariate comparison
### C) "proj_l48" an 'sf' object of the 2018 conterminous U.S.
source("code/preparation.R") # Note: it will take many minutes to run the Lee's L statistics

######################
# SUBPLOT GENERATION #
######################

# graphical expansion factor (to easily increase png resolution)
f <- 0.33

# Lee's L for Male Lung Cancer Mortality Rate by Male Ever Smoking Prevalence
## Differences from manuscript figures: Smaller graphical parameters, no legend, no title, no axes
LeeHex <- ggplot2::ggplot() + 
  ggplot2::geom_sf(data = Lee,
                   ggplot2::aes(fill = LESM_LeeFDR),
                   color = "black",
                   lwd = 0.1*f) +
  ggplot2::geom_sf(data = proj_l48,
                   fill = "transparent",
                   color = "black",
                   lwd = 0.25*f) +
  ggplot2::scale_fill_manual(values = c("#b2182b", "#b29f18", "#18b29f", "#2166ac", "grey90", "white"),
                             drop = FALSE,
                             na.translate = FALSE) + 
  ggplot2::guides(fill = "none") + 
  ggplot2::theme_void()

#####################
# CREATE HEXSTICKER #
#####################

# Additional Packages
library(hexSticker)

s <- hexSticker::sticker(subplot = LeeHex,
                         package = "Geographic\nPatterns in\nU.S. Lung Cancer\nMortality and\nCigarette Smoking",
                         p_size = 2.5, p_x = 1.37, p_y = 0.95, p_color = "black", # title
                         s_x = 0.05, s_y = 0.982, s_width = 2.36, s_height = 2.36, # symbol
                         h_fill = "white", # inside
                         h_color = "black", # outline
                         dpi = 1000, # resolution
                         filename = "hex/hex.png",
                         white_around_sticker = F)

# -------------------------------- END OF CODE --------------------------------- #
