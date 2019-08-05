# genscore
This repository contains the generalized score matching estimator introduced in the paper "Generalized Score Matching for Non-Negative Data" (http://www.jmlr.org/papers/volume20/18-278/18-278.pdf), an estimator for high-dimensional graphical models or parameters in truncated distributions. It is a generalization of the regularized score matching estimator in "Estimation of High-Dimensional Graphical Models Using Regularized Score Matching" (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5476334/, https://github.com/linlina/highscore).

The code here covers 
1. the univariate truncated normal distribution, 
1. Gaussian graphical models, 
1. truncated Gaussian graphical models, 
1. exponential square-root graphical models (Inouye et al, 2016),
1. "gamma graphical models" (Yu et al, 2019),
1. "a-b models" (Yu et al, 2019).

# Installation
install.packages("devtools")
devtools::install_github("sqyu/genscore")

# References
Some parts of the code were initially dervied from https://github.com/linlina/highscore.
