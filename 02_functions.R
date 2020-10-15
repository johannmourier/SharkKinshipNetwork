## Customized functions used in the paper: Specialized foraging underlies dolphin active social preferences

# Ordering matrix cf. SOCPROG ----
order_matrix <- function(x){
  newMATRIX <- x[sort(rownames(x), decreasing = FALSE),
                 sort(colnames(x), decreasing = FALSE)]
}

# Unfold Matrix ----
matrix_unfold <- function(x) {
  x[lower.tri(x, diag=FALSE)]
}

# Build Matrix from Unfolded matrix ----
matrix_folding <- function(matrix.unfold, matrix.labels, matrix.reference){
  tmp.matrix <- matrix(0, 
                       nrow = dim(matrix.reference)[1], 
                       ncol = dim(matrix.reference)[2],
                       dimnames = list(matrix.labels, matrix.labels)) 
  tmp.matrix[lower.tri(tmp.matrix)] <- matrix.unfold 
  return(tmp.matrix)
}

# SRI (by Mauricio Cantor) ----
SRI =  function (matr) {
  if (any(is.na(matr))) {
    matr <- na.omit(matr)
    cat("The data matrix contains NA, and have been removed.\n")
  }
  
  matr1 = matr
  N <- nrow(matr1)
  matr1[matr1 > 1] <- 1
  n <- apply(matr1, 2, sum)
  tmatr <- t(matr1)
  df <- as.matrix(t(matr))
  a <- df %*% t(df)
  b <- df %*% (1 - t(df))
  c <- (1 - df) %*% t(df)
  d <- ncol(df) - a - b - c
  
  Dice <- data.frame()
  inmat <- data.frame()
  denmat <- data.frame()
  
  for (i in 1:nrow(a)) {
    for (j in 1:ncol(a)) {
      # Simple Ratio Index
      Dice[i, j] <- 1 * a[i, j]/(1 * a[i, j] + b[i, j] + c[i, j])
      
      # Numerator of the SRI
      inmat[i, j] <- 1 * a[i, j] 
      
      # Denominator of the SRI
      denmat[i, j] <- (1 * a[i, j] + b[i, j] + c[i, j])
    }
  }
  
  rownames(Dice)=colnames(Dice)=colnames(matr)
  rownames(inmat)=colnames(inmat)=colnames(matr)
  rownames(denmat)=colnames(denmat)=colnames(matr)
  
  # Returns a list of 3 levels
  list(SRI = Dice, # Simple-Ratio indices
       SRI.numerator = inmat,  # Numerator of the SRI
       SRI.denominator = denmat) # Denominator of the SRI
}

# Permute GBI Matrix (by Mauricio Cantor) ----
null_checkerboard <- function(mat.gbi, iter, ...) {
  tmp.permute <- permatswap(mat.gbi, 
                            times = iter, 
                            method = "quasiswap", 
                            fixedmar = "both", 
                            shuffle = "both", 
                            mtype = "prab")
}


# Subset Matrix ----
matrix_subset <- function(x, info_subset) {
  x[which(rownames(x) %in% info_subset),
    which(colnames(x) %in% info_subset)]
}

# get giant component of a graph (by Mauricio Cantor) ----
giant.component <- function(graph) { 
  cl <- clusters(graph) 
  induced.subgraph(graph, which(cl$membership == which.max(cl$csize)))
} 

# Make matrices symmetric
mat_sym <- function(matrix) {
  matrix[upper.tri(matrix)] <- t(matrix)[lower.tri(matrix)]
  return(matrix)
}

# Change the number os breaks in ggplot2
# https://stackoverflow.com/questions/28436855/change-the-number-of-breaks-using-facet-grid-in-ggplot2

equal_breaks <- function(n = 3, s = 0.05, ...){
  function(x){
    # rescaling
    d <- s * diff(range(x)) / (1+2*s)
    seq(min(x)+d, max(x)-d, length=n)
  }
}

# Labels in 10^ format
# https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

# Rescale
#' @title Rescaling distribution
#' @description Rescales a distribution to a chosen range
#' @param x vector to be rescaled
#' @param r.out vector with limits of the new distribution
#' @return vector with rescaled distribution
#' @examples
#' # generating 10-sample uniform distribution from 1 to 100 and rescaling it to 0 to 1
#' data <- runif(10, 1, 100)
#' rescale(data, r.out=c(0,1))

rescale <- function(x, r.out) {
  p <- (x - min(x)) / (max(x) - min(x))
  r.out[[1]] + p * (r.out[[2]] - r.out[[1]])
}


# Theme for ggplot2 ----
theme_violin <- theme(strip.background = element_blank(), 
                      strip.placement = "outside",
                      strip.text = element_blank(),
                      axis.text = element_text(size = 10, colour = "black"),
                      axis.title = element_text(size = 10, colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.line.x = element_line(colour = "black"),
                      axis.title.x=element_blank(),
                      legend.position = "none",
                      panel.spacing = unit(1, "lines"),
                      plot.margin = unit(c(0, 0, 0.2, 0.1), "cm"))