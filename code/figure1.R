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
# A) Code to generate Figure 1 in the manuscript
# B) Male Lung Cancer Mortality Rate by Male Ever Smoking Prevalence
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

# Male Lung Cancer Mortality Rate
## Preparation
LCM <- Lee[stats::complete.cases(Lee$Male_Rate), ]
bks <- as.vector(quantile(LCM$Male_Rate, probs = seq(0, 1, 0.1), na.rm = TRUE))
bks[1] <- 0 # set lowest limit to zero
LCM$at <- cut(LCM$Male_Rate, bks, dig.lab = 20)
levels(LCM$at) <- levels(LCM$at)[gtools::mixedorder(levels(LCM$at), decreasing = TRUE)]
levels(LCM$at) <- stringr::str_replace(levels(LCM$at), ",", " - ") %>%
  stringr::str_replace(., "[(]", "") %>%
  stringr::str_replace(., "[]]", "")
tmp <- sprintf("%.1f", as.numeric(unlist(strsplit(levels(LCM$at), "-"))))
levels(LCM$at) <- paste(tmp[seq(1, 19, 2)], tmp[seq(2, 20, 2)], sep = " - ")
levels(LCM$at) <- c(levels(LCM$at), "Suppressed")
levels(LCM$at) <- levels(LCM$at)[gtools::mixedorder(levels(LCM$at), decreasing = FALSE)]
LCM$at <- relevel(LCM$at, "Suppressed")
levels(LCM$at) <- c("Suppressed", "13.4 - 58.8", "58.9 - 67.8", "67.9 - 76.0", "76.1 - 83.0", "83.1 - 90.0",
                           "90.1 - 97.2", "97.3 - 105.1", "105.2 - 114.3", "114.4 - 126.8", "126.9 - 297.3")
## Plotting
CMRM <- ggplot2::ggplot() + 
  ggplot2::geom_sf(data = Lee,
                   fill = "white",
                   color = "black") +
  ggplot2::geom_sf(data = LCM,
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

# Male Ever Smoking Prevalence
## Preparation
ESM <- Lee[stats::complete.cases(Lee$ESM_Pct_1997_2003), ]
bks <- as.vector(quantile(ESM$ESM_Pct_1997_2003, probs = seq(0, 1, 0.1), na.rm = TRUE))
bks[1] <- 0 # set lowest limit to zero
ESM$at <- cut(ESM$ESM_Pct_1997_2003, bks, dig.lab = 20)
levels(ESM$at) <- levels(ESM$at)[gtools::mixedorder(levels(ESM$at), decreasing = TRUE)]
levels(ESM$at) <- stringr::str_replace(levels(ESM$at), ",", " - ") %>%
  stringr::str_replace(., "[(]", "") %>%
  stringr::str_replace(., "[]]", "")
tmp <- sprintf("%.1f", as.numeric(unlist(strsplit(levels(ESM$at), "-"))))
levels(ESM$at) <- paste(tmp[seq(1, 19, 2)], tmp[seq(2, 20, 2)], sep = " - ")
levels(ESM$at) <- c(levels(ESM$at), "Suppressed")
levels(ESM$at) <- levels(ESM$at)[gtools::mixedorder(levels(ESM$at), decreasing = FALSE)]
ESM$at <- relevel(ESM$at, "Suppressed")
levels(ESM$at) <- c("Suppressed", "21.6 - 50.7", "50.8 - 54.0", "54.1 - 56.1", "56.2 - 57.7", "57.8 - 59.2",
                           "59.3 - 60.5", "60.6 - 61.8", "61.9 - 63.4", "63.5 - 65.2", "65.3 - 74.6")
## Plotting
PESM <- ggplot2::ggplot() + 
  ggplot2::geom_sf(data = Lee,
                   fill = "white",
                   color = "black") +
  ggplot2::geom_sf(data = ESM,
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

# Lee's L for Male Lung Cancer Mortality Rate by Male Ever Smoking Prevalence
LESM <- ggplot2::ggplot() + 
  ggplot2::geom_sf(data = Lee,
                   ggplot2::aes(fill = LESM_LeeFDR),
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
# FIGURE 1 #
############

# Combine for Figure 1
BESM <- cowplot::align_plots(CMRM, PESM, LESM, align = "hv", axis = "tblr")
gg_inset_map <- cowplot::ggdraw() +
  cowplot::draw_plot(BESM[[1]], y =  0.33, scale = 0.9) +
  cowplot::draw_plot(BESM[[2]], y =  0.00, scale = 0.9) +
  cowplot::draw_plot(BESM[[3]], y = -0.33, scale = 0.9)

# Print
grDevices::png(file = "figures/figure1.png", width = 850*f, height = 1100*f)
print(gg_inset_map)
grDevices::dev.off()

# ----------------------------- END OF CODE ----------------------------- #
