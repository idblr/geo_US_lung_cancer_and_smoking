# ----------------------------------------------------------------------- #
# Geographic Patterns in U.S. Lung Cancer Mortality and Cigarette Smoking
# ----------------------------------------------------------------------- #
#
# Created by: Ian Buller, Ph.D., M.A. (GitHub: @idblr)
# Created on: January 15, 2022
#
# Most recently modified by:
# Most recently modified on:
#
# Notes:
# A) Code to generate Supplemental Figure 2 in the manuscript
# B) Female Lung Cancer Mortality Rate by Female Current Smoking Prevalence
# C) For consistency with the manuscript, must run all four Lee's L statistics in order with the RNG seed (see 'preparation.R' file)
# ----------------------------------------------------------------------- #

####################
# DATA IMPORTATION #
####################

# Use the code found in 'functions.R' and 'preparation.R' files
## Loads three objects
### A) 'Lee' an 'sf' object of the U.S. county-level lung cancer mortality rates, smoking prevalences, and Lee's L values
### B) 'FDRpvals' a 'data.frame' object of significant p-value cutoff for each bivariate comparison
### C) "proj_l48_coast" an 'sf' object of the 2018 conterminous U.S.
source("code/functions.R")
source("code/preparation.R") # Note: it will take many minutes to run the Lee's L statistics

############
# SETTINGS #
############

# graphical expansion factor (to easily increase png resolution)
f <- 5 

## Colors for decile colorkey
brks <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")

#########################
# MULTI-PLOT COMPONENTS #
#########################

# Female Lung Cancer Mortality Rate
## Preparation
LCF <- Lee[stats::complete.cases(Lee$Female_Rate), ]
bks <- as.vector(quantile(LCF$Female_Rate, probs = seq(0, 1, 0.1), na.rm = TRUE))
bks[1] <- 0 # set lowest limit to zero
LCF$at <- cut(LCF$Female_Rate, bks, dig.lab = 20)
levels(LCF$at) <- levels(LCF$at)[gtools::mixedorder(levels(LCF$at), decreasing = TRUE)]
levels(LCF$at) <- stringr::str_replace(levels(LCF$at), ",", " - ") %>%
  stringr::str_replace(., "[(]", "") %>%
  stringr::str_replace(., "[]]", "")
tmp <- sprintf("%.1f", as.numeric(unlist(strsplit(levels(LCF$at), "-"))))
levels(LCF$at) <- paste(tmp[seq(1, 19, 2)], tmp[seq(2, 20, 2)], sep = " - ")
levels(LCF$at) <- c(levels(LCF$at), "Suppressed")
levels(LCF$at) <- levels(LCF$at)[gtools::mixedorder(levels(LCF$at), decreasing = FALSE)]
LCF$at <- relevel(LCF$at, "Suppressed")
levels(LCF$at) <- c("Suppressed", "12.4 - 38.7", "38.8 - 44.1", "44.2 - 48.0", "48.1 - 51.6", "51.7 - 55.1",
                           "55.2 - 58.3", "58.4 - 61.8", "61.9 - 66.0", "66.1 - 72.5", "72.6 - 133.7")
## Plotting
CMRF <- ggplot2::ggplot() + 
  ggplot2::geom_sf(data = Lee,
                   fill = "white",
                   color = "black") +
  ggplot2::geom_sf(data = LCF,
                   ggplot2::aes(fill = at),
                   color = "black") +
  ggplot2::geom_sf(data = proj_l48_coast,
                   fill = "transparent",
                   color = "black",
                   lwd = 0.33*f) +
  ggplot2::scale_fill_manual(values = c("white", brks),
                             drop = FALSE) + 
  ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5*f),
                                               title = "Mortality Rate",
                                               order = 1,
                                               reverse = TRUE)) +
  ggplot2::labs(title = "(A)") +
  ggplot2::theme_minimal(base_size = 14*f)

# Female Current Smoking Prevalence
## Preparation
CSF <- Lee[stats::complete.cases(Lee$CSF_Pct_1997_2003), ]
bks <- as.vector(quantile(Lee$CSF_Pct_1997_2003, probs = seq(0, 1, 0.1), na.rm = TRUE))
bks[1] <- 0 # set lowest limit to zero
CSF$at <- cut(CSF$CSF_Pct_1997_2003, bks, dig.lab = 20)
levels(CSF$at) <- levels(CSF$at)[gtools::mixedorder(levels(CSF$at), decreasing = TRUE)]
levels(CSF$at) <- stringr::str_replace(levels(CSF$at), ",", " - ") %>%
  stringr::str_replace(., "[(]", "") %>%
  stringr::str_replace(., "[]]", "")
tmp <- sprintf("%.1f", as.numeric(unlist(strsplit(levels(CSF$at), "-"))))
levels(CSF$at) <- paste(tmp[seq(1, 19, 2)], tmp[seq(2, 20, 2)], sep = " - ")
levels(CSF$at) <- c(levels(CSF$at), "Suppressed")
levels(CSF$at) <- levels(CSF$at)[gtools::mixedorder(levels(CSF$at), decreasing = FALSE)]
CSF$at <- relevel(CSF$at, "Suppressed")
levels(CSF$at) <- c("Suppressed", "1.7 - 17.0", "17.1 - 19.4", "19.5 - 21.0", "21.1 - 22.2", "22.3 - 23.3",
                           "23.4 - 24.4", "24.5 - 25.6", "25.7 - 27.0", "27.1 - 29.1", "29.2 - 48.4")
## Plotting
PCSF <- ggplot2::ggplot() + 
  ggplot2::geom_sf(data = Lee,
                   fill = "white",
                   color = "black") +
  ggplot2::geom_sf(data = CSF,
                   ggplot2::aes(fill = at),
                   color = "black") +
  ggplot2::geom_sf(data = proj_l48_coast,
                   fill = "transparent",
                   color = "black",
                   lwd = 0.33*f) +
  ggplot2::scale_fill_manual(values = c("white", brks),
                             drop = FALSE) + 
  ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5*f),
                                               title = "Smoking Prevalence",
                                               order = 1,
                                               reverse = TRUE)) + 
  ggplot2::labs(title = "(B)") +
  ggplot2::theme_minimal(base_size = 14*f)

# Lee's L for Female Lung Cancer Mortality Rate by Female Current Smoking Prevalence
LCSF <- ggplot2::ggplot() + 
  ggplot2::geom_sf(data = Lee,
                   ggplot2::aes(fill = LCSF_LeeFDR),
                   color = "black") +
  ggplot2::geom_sf(data = proj_l48_coast,
                   fill = "transparent",
                   color = "black",
                   lwd = 0.33*f) +
  ggplot2::scale_fill_manual(values = c("#b2182b", "#b29f18", "#18b29f", "#2166ac", "grey90", "white"),
                             drop = FALSE,
                             na.translate = FALSE) +
  ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5*f),
                                               title = "Lee's L Clusters")) + 
  ggplot2::labs(title = "(C)") +
  ggplot2::theme_minimal(base_size = 14*f)

#########################
# SUPPLEMENTAL FIGURE 2 #
#########################

# Combine for Supplemental Figure 2
BCSF <- cowplot::align_plots(CMRF, PCSF, LCSF, align = "hv", axis = "tblr")
gg_inset_map <- cowplot::ggdraw() +
  cowplot::draw_plot(BCSF[[1]], y =  0.33, scale = 0.9) +
  cowplot::draw_plot(BCSF[[2]], y =  0.00, scale = 0.9) +
  cowplot::draw_plot(BCSF[[3]], y = -0.33, scale = 0.9)

# Print
grDevices::png(file = "figures/supplemental2.png", width = 850*f, height = 1100*f)
print(gg_inset_map)
grDevices::dev.off()

# ----------------------------- END OF CODE ----------------------------- #
