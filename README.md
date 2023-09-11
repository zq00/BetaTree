# BetaTree

## Overview

This R package constructs a Beta tree histogram for multivariate data by growing a k-d tree on the order statistics, computing simultaneous confidence bound of average density in each region, and selecting the largest regions where data are approximately uniform. The Beta tree histogram provides a succinct description of the data, as well as simultaneous confidence intervals of the average density in each region. 

## Getting started

- Package website: https://zq00.github.io/BetaTrees/ 
- You can find the package **vignettes** under the Articles tab. 
- You can read more about the Beta Tree histogram in the paper.  

## Installation

You can install the package using 

```R
install.packages("devtools")
devtools::install_github("zq00/BetaTrees")
```

## Function documentation

You can find the function documentations under the  Reference tab. To get started, you can take a look at the function `BuilgHist()`, which computes a Beta tree histogram.

## Source code

The source code is located at the [Github -> R](https://github.com/zq00/BetaTrees/tree/main/R) folder. 

## Feedback

If you encounter error or would like to provide feedback, please use [Github -> Issues](https://github.com/zq00/BetaTrees/issues) to reach us. Thank you! 
