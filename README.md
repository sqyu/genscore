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

# Installation from CRAN
install.packages("genscore")

# Installation from GitHub
```R
install.packages(c("devtools", "knitr"))
devtools::install_github("sqyu/genscore", build_vignettes=TRUE)
# Set build_vignettes to FALSE if you do not wish to build the vignette (which takes a few minutes).
```

# Usage
For a complete guide to its usage, please consult the vignette [here](vignettes/gen_vignette.Rmd) (or [here](https://htmlpreview.github.io/?https://github.com/sqyu/genscore/blob/master/vignettes/gen_vignette.html) for the precompiled html).
```R
vignette("gen_vignette")
```

# References
Some parts of the code were initially dervied from https://github.com/linlina/highscore and http://www1.maths.leeds.ac.uk/~wally.gilks/adaptive.rejection/web_page/Welcome.html.

John Aitchison. A general class of distributions on the simplex. Journal of the Royal Statistical Society: Series B (Methodological), 47(1):136–146, 1985. <https://doi.org/10.1111/j.2517-6161.1985.tb01341.x>

David Inouye, Pradeep Ravikumar, and Inderjit Dhillon. Square root graphical models: Multivariate generalizations of univariate exponential families that permit positive de- pendencies. In Proceedings of The 33rd International Conference on Machine Learning, volume 48 of Proceedings of Machine Learning Research, pages 2445–2453, 2016. <http://proceedings.mlr.press/v48/inouye16.html>

Lina Lin, Mathias Drton, and Ali Shojaie. Estimation of high-dimensional graphical models using regularized score m    atching. Electron. J. Stat., 10(1):806–854, 2016. <https://doi.org/10.1214/16-EJS1126>

Shiqing Yu, Mathias Drton, and Ali Shojaie. Graphical models for non-negative data using generalized score matching. In International Conference on Artificial Intelligence and Statistics, pages 1781–1790, 2018. <http://proceedings.mlr.press/v84/yu18b.html>

Shiqing Yu, Mathias Drton, and Ali Shojaie. Generalized score matching for non-negative data. Journal of Machine Learning Research, 20(76):1–70, 2019. <http://jmlr.org/papers/v20/18-278.html>





