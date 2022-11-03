## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = F,
  message = F
)

## ----setup, include = F-------------------------------------------------------
library(tidyverse)
devtools::load_all(path = "~/Desktop/Topics/kdhist/package/BetaTrees")

## -----------------------------------------------------------------------------
n <- 2000 # num of obs.
d <- 2 # num of dim.

## ---- echo = F----------------------------------------------------------------
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

## ---- echo = F, fig.height = 3, fig.width= 4, fig.align='center'--------------
ggplot() + 
  geom_point(aes(x = X[,1], y = X[,2] , color = as.factor(N)), size = 0.5) + 
  theme_bw() + 
  labs(color = "Component", x = "X", y="Y") + 
  theme(legend.position = c(0.85, 0.8)) 

## -----------------------------------------------------------------------------
hist <- BuildHist(X, alpha = 0.1, method = "weighted_bonferroni")

## -----------------------------------------------------------------------------
hist[1,]

## ---- fig.height = 3, fig.width= 5, fig.align='center'------------------------
hist <- BuildHist(X, alpha = 0.1, method = "weighted_bonferroni", plot = TRUE)

## ---- fig.height = 3, fig.width= 5, fig.align='center'------------------------
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
ahat <- ConfLevel(tree$nd, alpha = 0.01, method = "weighted_bonferroni")
ahat

## -----------------------------------------------------------------------------
tree <- SetBounds(tree$kdtree, ahat, n)

## -----------------------------------------------------------------------------
tree$leftchild$leftchild$rightchild$rightchild$leftchild$leftchild

## -----------------------------------------------------------------------------
B <- matrix(nrow = 0, ncol = (2*d + 5)) 
hist <- SelectNodes(tree, B, ahat, n)

