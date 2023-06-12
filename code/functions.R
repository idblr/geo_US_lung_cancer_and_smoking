# ----------------------------------------------------------------------- #
# Geographic Patterns in U.S. Lung Cancer Mortality and Cigarette Smoking
# ----------------------------------------------------------------------- #
#
# Created by: Ian Buller, Ph.D., M.A. (GitHub: @idblr)
# Created on: January 15, 2022
#
# Most recently modified by: @idblr
# Most recently modified on: June 5, 2023
#
# Notes:
# A) Code for essential functions to calculate the Lee's L statistic and False Discovery Rate
# B) 05/17/2022 (@idblr): Consistent variable name based on merged data object; rename NAs as "Suppressed" in LeeL() function
# C) 05/17/2022 (@idblr): Correctly set "dat" in LeeL() function
# D) 06/05/2023 (@idblr): Update to FDR calculation for multiple testing correction
# ----------------------------------------------------------------------- #

cat("\nLoading custom functions for Lee's L statistic\n")

# Select not in
`%notin%` <- Negate(`%in%`) # https://www.r-bloggers.com/the-notin-operator/

# False Discovery Rate (Benjamini & Hochberg, 1995; DOI: 10.1111/j.2517-6161.1995.tb02031.x)
fdr <- function(pvals, alpha) {
  pcrit <- NULL
  m <- length(pvals)
  for (i in 1:m) {
    if (pvals[i] <= (i/m) * alpha) { pcrit <- pvals[i] }
  }
  return(max(pcrit, pvals[1]))
}

# Local Lee
## Permutations for the Bivariate Moran's I
simula_lee <- function(x, y, listw, nsim = nsim, zero.policy = NULL, na.action = na.fail){
  
  if (deparse(substitute(na.action)) == "na.pass") 
    stop("na.pass not permitted")
  na.act <- attr(na.action(cbind(x, y)), "na.action")
  x[na.act] <- NA
  y[na.act] <- NA
  x <- na.action(x)
  y <- na.action(y)
  if (!is.null(na.act)) {
    subset <- !(1:length(listw$neighbours) %in% na.act)
    listw <- subset(listw, subset, zero.policy = zero.policy)
  }
  n <- length(listw$neighbours)
  if ((n != length(x)) | (n != length(y))) 
    stop("objects of different length")
  gamres <- suppressWarnings(nsim > gamma(n + 1))
  if (gamres) 
    stop("nsim too large for this number of observations")
  if (nsim < 1) 
    stop("nsim too small")
  xy <- data.frame(x, y)
  S2 <- sum((unlist(lapply(listw$weights, sum)))^2)
  
  lee_boot <- function(var, i, ...) {
    return(spdep::lee(x = var[i, 1], y = var[i, 2], ...)$localL)
  }
  
  res <- boot::boot(xy, statistic = lee_boot, R = nsim, sim = "permutation", 
                    listw = listw, n = n, S2 = S2, zero.policy = zero.policy)
}

## Calculate Local Lee's L statistic with multiple testing correction by False Discovery Rate
LeeL <- function(dat, x, y, label, numsim = 100,...) {
  
  if(class(dat) != "SpatialPolygonsDataFrame") stop("dat is not of class 'SpatialPolygonsDataFrame'")
  
  temp <- dat[stats::complete.cases(dat[[y]]) & stats::complete.cases(dat[[x]]), ]
  X <- temp[[x]]
  Y <- temp[[y]]
  
  nb0 <- spdep::poly2nb(temp) # QUEEN
  no_neighs <- which(spdep::card(nb0) == 0) # empty geometries
  
  ### Nearest neighbor for empty geometries
  #### https://stat.ethz.ch/pipermail/r-sig-geo/2013-July/018864.html
  k1nb <- spdep::knn2nb(spdep::knearneigh(sp::coordinates(temp),
                                          k = 1,
                                          longlat = !sp::is.projected(temp)))
  res <- k1nb[no_neighs]
  nb0[no_neighs] <- res
  attr(nb0, "sym") <- spdep::is.symmetric.nb(nb0, force = TRUE)
  #### re-assign the now incorrect symmetry attribute
  nb0 <- spdep::make.sym.nb(nb0)
  
  ### Adjacency Matrix (Queen)
  col.W <- spdep::nb2listw(nb0, style = "W")
  lw <- spdep::nb2listw(nb0, style = "B", zero.policy = T)
  W  <- as(lw, "symmetricMatrix")
  W  <- as.matrix(W / Matrix::rowSums(W))
  W[which(is.na(W))] <- 0
  
  ### Calculating the index and its simulated distribution for local values
  m <- spdep::lee(x = X,
                  y = Y,
                  listw = col.W,
                  n = length(X),
                  zero.policy = TRUE,
                  NAOK = TRUE)
  
  local_sims <- simula_lee(x = X,
                           y = Y,
                           listw = col.W,
                           nsim = numsim,
                           zero.policy = TRUE,
                           na.action = na.omit)
  
  m_i <- m[[2]]  # local values
  
  ### Identifying the significant values 
  #### significant based on empirical p-value
  pvals <- (rowSums(t(abs(local_sims$t)) >= abs(local_sims$t0))+1)/((nrow(local_sims[[2]]))+1)
  
  #### Correction for Multiple Testing
  ##### False Discovery Rate
  sort_pvals <- sort(as.vector(pvals))
  out_alpha <- fdr(sort_pvals, alpha)
  sig <- pvals < out_alpha
  
  ### Preparing for plotting
  temp$sig <- sig
  temp$pval <- pvals
  
  #### Identifying the Lee's L patterns
  Xp <- scale(X)
  Yp <- scale(Y)
  patterns <- as.character( interaction(Xp > 0, W %*% Yp > 0) ) 
  patterns <- patterns |>
    stringr::str_replace_all("TRUE", "High") |>
    stringr::str_replace_all("FALSE", "Low")
  patterns[temp$sig == 0] <- "Not significant"
  temp$patterns <- patterns
  temp$out <- factor(temp$patterns,
                      levels = c("High.High", "Low.High", "High.Low",
                                 "Low.Low", "Not significant", "NA"),
                      labels = c("High smoking - High mortality", "Low smoking - High mortality",
                                 "High smoking - Low mortality", "Low smoking - Low mortality",
                                 "Not significant", "Suppressed"))
  
  tmp <- temp[, c("GEOID", "STATEFP", "COUNTYFP", "NAME", "pval", "out")]
  names(tmp) <- c("GEOID", "STATEFP", "COUNTYFP", "NAME", paste(label, "_LeeP", sep = ""), paste(label, "_LeeFDR", sep = ""))
  out <- merge(dat, tmp)
  return(list(out = out,
              FDRpval = out_alpha))
}

# ----------------------------- END OF CODE ----------------------------- #
