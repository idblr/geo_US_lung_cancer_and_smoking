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
# A) Code for essential functions to calculate the Lee's L statistic and False Discovery Rate
# ----------------------------------------------------------------------- #

cat("\nLoading custom functions for Lee's L statistic\n")

# Select not in
`%notin%` <- Negate(`%in%`) # https://www.r-bloggers.com/the-notin-operator/

# False Discovery Rate
fdr <- function(pvals, alpha = 0.05) {
  m <- length(pvals)
  for (i in 1:length(pvals)) {
    if (pvals[i] <= (i/m) * alpha) { return(pvals[i]) } 
  }
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
  
  if(class(alcove_proj) != "SpatialPolygonsDataFrame") stop("dat is not of class 'SpatialPolygonsDataFrame'")
  
  temp <- dat[stats::complete.cases(dat[[y]]) & stats::complete.cases(dat[[x]]), ]
  X <- temp[[x]]
  Y <- temp[[y]]
  
  nb0 <- spdep::poly2nb(temp) # QUEEN
  no_neighs <- which(spdep::card(nb0) == 0) # 3 empty geometries
  
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
  sort_pvals <- sort(as.vector(pvals), decreasing = TRUE)
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
  
  tmp <- temp[, c("GEOID.1", "STATEFP", "COUNTYFP", "NAME.1", "pval", "out")]
  names(tmp) <- c("GEOID.1", "STATEFP", "COUNTYFP", "NAME", paste(label, "_LeeP", sep = ""), paste(label, "_LeeFDR", sep = ""))
  out <- merge(dat, tmp)
  return(list(out = out,
              FDRpval = out_alpha))
}

# ----------------------------- END OF CODE ----------------------------- #
