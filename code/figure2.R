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
# A) Code to generate Figure 2 in the manuscript
# B) Female Lung Cancer Mortality Rate by Female Ever Smoking Prevalence
# C) For consistency with the manuscript, must run all four Lee's L statistics in order with the RNG seed (see 'preparation.R' file)
# D) 05/17/2022: Update for correctly suppressed counties for each result by sex
# ----------------------------------------------------------------------- #

####################
# DATA IMPORTATION #
####################

# Use the code found in 'functions.R' and 'preparation.R' files
## Loads three objects
### A) 'Lee' an 'sf' object of the U.S. county-level lung cancer mortality rates, smoking prevalences, and Lee's L values
### B) 'FDRpvals' a 'data.frame' object of significant p-value cutoff for each bivariate comparison
### C) "proj_l48" an 'sf' object of the 2018 conterminous U.S.
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
ESF <- Lee[stats::complete.cases(Lee$Female_Rate), ]
bks <- as.vector(quantile(ESF$Female_Rate, probs = seq(0, 1, 0.1), na.rm = TRUE))
bks[1] <- 0 # set lowest limit to zero
ESF$at <- cut(ESF$Female_Rate, bks, dig.lab = 20)
levels(ESF$at) <- levels(ESF$at)[gtools::mixedorder(levels(ESF$at), decreasing = TRUE)]
levels(ESF$at) <- stringr::str_replace(levels(ESF$at), ",", " - ") %>%
  stringr::str_replace(., "[(]", "") %>%
  stringr::str_replace(., "[]]", "")
tmp <- sprintf("%.1f", as.numeric(unlist(strsplit(levels(ESF$at), "-"))))
levels(ESF$at) <- paste(tmp[seq(1, 19, 2)], tmp[seq(2, 20, 2)], sep = " - ")
levels(ESF$at) <- c(levels(ESF$at), "Suppressed")
levels(ESF$at) <- levels(ESF$at)[gtools::mixedorder(levels(ESF$at), decreasing = FALSE)]
ESF$at <- relevel(ESF$at, "Suppressed")
levels(ESF$at) <- c("Suppressed", "12.4 - 38.7", "38.8 - 44.1", "44.2 - 48.0", "48.1 - 51.6", "51.7 - 55.1",
                           "55.2 - 58.3", "58.4 - 61.8", "61.9 - 66.0", "66.1 - 72.5", "72.6 - 133.7")
## Plotting
CMRF <- ggplot2::ggplot() + 
  ggplot2::geom_sf(data = Lee,
                   fill = "white",
                   color = "black") +
  ggplot2::geom_sf(data = ESF,
                   ggplot2::aes(fill = at),
                   color = "black") +
  ggplot2::geom_sf(data = proj_l48,
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

# Female Ever Smoking Prevalence
## Preparation
ESF <- Lee[stats::complete.cases(Lee$ESF_Pct_1997_2003), ]
bks <- as.vector(quantile(ESF$ESF_Pct_1997_2003, probs = seq(0, 1, 0.1), na.rm = TRUE))
bks[1] <- 0 # set lowest limit to zero
ESF$at <- cut(ESF$ESF_Pct_1997_2003, bks, dig.lab = 20)
levels(ESF$at) <- levels(ESF$at)[gtools::mixedorder(levels(ESF$at), decreasing = TRUE)]
levels(ESF$at) <- stringr::str_replace(levels(ESF$at), ",", " - ") %>%
  stringr::str_replace(., "[(]", "") %>%
  stringr::str_replace(., "[]]", "")
tmp <- sprintf("%.1f", as.numeric(unlist(strsplit(levels(ESF$at), "-"))))
levels(ESF$at) <- paste(tmp[seq(1, 19, 2)], tmp[seq(2, 20, 2)], sep = " - ")
levels(ESF$at) <- c(levels(ESF$at), "Suppressed")
levels(ESF$at) <- levels(ESF$at)[gtools::mixedorder(levels(ESF$at), decreasing = FALSE)]
ESF$at <- relevel(ESF$at, "Suppressed")
levels(ESF$at) <- c("Suppressed", "1.5 - 33.4", "33.5 - 36.7", "36.8 - 38.8", "38.9 - 40.5", "40.6 - 42.1",
                           "42.2 - 43.7", "43.8 - 45.2", "45.3 - 46.7", "46.8 - 49.0", "49.1 - 60.4")
## Plotting
PESF <- ggplot2::ggplot() + 
  ggplot2::geom_sf(data = Lee,
                   fill = "white",
                   color = "black") +
  ggplot2::geom_sf(data = ESF,
                   ggplot2::aes(fill = at),
                   color = "black") +
  ggplot2::geom_sf(data = proj_l48,
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

# Lee's L for Female Lung Cancer Mortality Rate by Female Ever Smoking Prevalence
LESF <- ggplot2::ggplot() + 
  ggplot2::geom_sf(data = Lee,
                   ggplot2::aes(fill = LESF_LeeFDR),
                   color = "black") +
  ggplot2::geom_sf(data = proj_l48,
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

############
# FIGURE 2 #
############

# Combine for Figure 2
BESF <- cowplot::align_plots(CMRF, PESF, LESF, align = "hv", axis = "tblr")
gg_inset_map <- cowplot::ggdraw() +
  cowplot::draw_plot(BESF[[1]], y =  0.33, scale = 0.9) +
  cowplot::draw_plot(BESF[[2]], y =  0.00, scale = 0.9) +
  cowplot::draw_plot(BESF[[3]], y = -0.33, scale = 0.9)

# Print
grDevices::png(file = "figures/figure2.png", width = 850*f, height = 1100*f)
print(gg_inset_map)
grDevices::dev.off()

# ----------------------------- END OF CODE ----------------------------- #
