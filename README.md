# pcglassoFast
`pcglassoFast`: Fast Partial Correlation Graphical LASSO

You can install it from a local directory:

```r
library(devtools)

install_local("/path/to/pcglassoFast")
```

Alternatively, if you have access to GitHub:

```r
library(devtools)

install_github("anonymous/pcglassoFast")
```

Note that you may need to install additional compilation tools to build the C++ and Fortran code included in the package. On Windows, this usually requires [Rtools](https://cran.r-project.org/bin/windows/Rtools/). On macOS, see the [R for macOS tools](https://mac.r-project.org/tools/) page.

## Example

```r
library(pcglassoFast)
library(MASS)

set.seed(1)

p <- 4
n <- 30
R.true <- toeplitz(c(c(1, -0.5), rep(0, p - 2)))
D.true <- sqrt(rchisq(p, 3))
K.true <- diag(D.true) %*% R.true %*% diag(D.true) # sparse precision matrix
S.true <- solve(K.true)

Z <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = S.true)
S <- cov(Z)

lambda <- 0.3
alpha <- 0
res <- pcglassoFast(S, lambda, alpha)

# Estimated precision matrix
res$Sinv

# True precision matrix
K.true
```

## **Report an issue**
Any bugs encountered when using the package can be reported [here](https://github.com/anonymous/pcglassoFast/issues).

## **Notes**
Part of the code (`ROptimDual()` and the Fortran code) was adapted from the `glassoFast` package: <https://github.com/JClavel/glassoFast>.

## Citation

If you use `pcglassoFast`, please cite both the package and the paper describing the method and algorithms:

```bibtex
@misc{anonymous2025pcglassofast,
  title  = {{pcglassoFast}: Fast Partial Correlation Graphical LASSO},
  author = {[Anonymous et al.]},
  year   = {2025},
  url    = {https://github.com/anonymous/pcglassoFast},
  note   = {R package}
}
```
