## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = F,
  message = F
)

## ----setup, include = F-------------------------------------------------------
library(tidyverse)
library(BetaTree)
library(mvtnorm)

## ----echo = F-----------------------------------------------------------------
n <- 2000 # num of obs.
d <- 2 # num of dim.

## ----echo = F-----------------------------------------------------------------
mu2d <- list(c(-1.5, 0.6), c(2, -1.5)) # mean 
sigma12d <- matrix(c(1, 0.5, 0.5, 1), d, d, byrow = T)
sigma22d <- matrix(c(1, 0, 0, 1), d, d, byrow = T)
sigma2d <- list(sigma12d, sigma22d) # covariance

N <- apply(rmultinom(n, 1, prob = c(0.4, 0.6)), 2, function(t) which(t == 1)) 
nn <- table(N) # number of elements from each mixture component
X <- matrix(0, n, d) # observations
# sample from each of the mixture
for(i in 1:2){
  R <- chol(sigma2d[[i]])
  G <- matrix(rnorm(nn[i] * 2, 0, 1), nn[i], d) %*% R
  X[which(N == i), ] <-   t(t(G) + mu2d[[i]])
}

## ----echo = F, fig.height = 3, fig.width= 4, fig.align='center'---------------
ggplot() + 
  geom_point(aes(x = X[,1], y = X[,2] , color = as.factor(N)), size = 0.5) + 
  theme_bw() + 
  labs(color = "Mixture", x = "X", y="Y") + 
  theme(legend.position = c(0.85, 0.8)) 

## ----fig.height = 3, fig.width= 5, fig.align='center'-------------------------
hist <- BuildHist(X, alpha = 0.1, method = "weighted_bonferroni", plot = T)

## -----------------------------------------------------------------------------
hist[1,]

## ----fig.height = 3, fig.width= 5, fig.align='center'-------------------------
hist <- BuildHist(X, alpha = 0.1, method = "weighted_bonferroni", bounded = T, 
                  option = "qt", q=c(0.05,0.05), 
                  plot = TRUE)

## -----------------------------------------------------------------------------
tree <- BuildKDTree(X, bounded = F)

## -----------------------------------------------------------------------------
tree$kdtree$leftchild$leftchild$rightchild$rightchild$leftchild$leftchild

## -----------------------------------------------------------------------------
tree$nd

## -----------------------------------------------------------------------------
ahat <- ConfLevel(tree$nd, alpha = 0.1, method = "weighted_bonferroni")
ahat

## -----------------------------------------------------------------------------
tree <- SetBounds(tree$kdtree, ahat, n)

## -----------------------------------------------------------------------------
tree$leftchild$leftchild$rightchild$rightchild$leftchild$leftchild

## -----------------------------------------------------------------------------
B <- matrix(nrow = 0, ncol = (2*d + 5)) 
hist <- SelectNodes(tree, B, ahat, n)

## ----simulate_data------------------------------------------------------------
n <- 10000
p <- 5
X <- matrix(rnorm(n*p), n, p)
X[,1:2] <- X[,1:2] %*% chol(matrix(c(1, 0.8, 0.8, 1), 2, 2))

## -----------------------------------------------------------------------------
alpha <- 0.1 
adaptive_beta_tree <- build_adaptive_histogram(X, 
                                               alpha  = alpha,
                                               thresh_marginal = alpha / ncol(X),
                                               thresh_interaction = qgamma(shape = ncol(X) - 1,rate = 1, p = 1 - alpha),
                                               method = "weighted_bonferroni")

## ----fig.height = 4, fig.width = 6--------------------------------------------
plot_dat <- plot_histogram_data( X = X,
                                 hist   = adaptive_beta_tree$hist,
                                 plot_coord = c(1, 2),
                                 nint  = 10,
                                 show_data  = TRUE,
                                 ndat     = 500)
plot_histogram(plot_dat)

