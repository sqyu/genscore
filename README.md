# genscore
This repository contains the generalized score matching estimator introduced in the paper "Generalized Score Matching for Non-Negative Data" (http://www.jmlr.org/papers/volume20/18-278/18-278.pdf), an estimator for high-dimensional graphical models or parameters in truncated distributions. It is a generalization of the regularized score matching estimator in "Estimation of High-Dimensional Graphical Models Using Regularized Score Matching" (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5476334/, https://github.com/linlina/highscore).

The current second version further generalizes the distributions to generalized domain types 
1. The real space,
1. The non-negative orthant of the real space,
1. A union of intervals as the uniform domain for each component,
1. The (p-1)-dimensional simplex (with all components positive and sum to 1), and
1. Intersections and unions of domains defined by polynomial inequalities.

The distributions covered include
1. the univariate truncated normal distribution, 
1. Gaussian graphical models, 
1. truncated Gaussian graphical models, 
1. exponential square-root graphical models (Inouye et al, 2016),
1. "gamma graphical models" (Yu et al, 2019),
1. "a-b models" (Yu et al, 2019), and
1. the A^d model (Aitchison, 1985).

# Installation
```R
install.packages(c("devtools", "knitr"))
devtools::install_github("sqyu/genscore", build_vignettes=TRUE)
```

# Usage
For a complete guide to its usage, please consult the vignette [here](vignettes/gen_vignette.Rmd) (or [here](https://htmlpreview.github.io/?https://github.com/sqyu/genscore/blob/master/vignettes/gen_vignette.html) for the precompiled html).
```R
vignette("gen_vignette")
```

# References
Some parts of the code were initially dervied from https://github.com/linlina/highscore.
