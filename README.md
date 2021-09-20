# DTAmc 

[![](https://www.r-pkg.org/badges/version/DTAmc?color=orange)](https://cran.r-project.org/package=DTAmc)
[![](https://img.shields.io/badge/devel%20version-0.0.8.9003-blue.svg)](https://github.com/maxwestphal/DTAmc)
[![R-CMD-check](https://github.com/maxwestphal/DTAmc/workflows/R-CMD-check/badge.svg)](https://github.com/maxwestphal/DTAmc/actions)
[![Codecov test coverage](https://codecov.io/gh/maxwestphal/DTAmc/branch/master/graph/badge.svg)](https://codecov.io/gh/maxwestphal/DTAmc?branch=master)
[![](https://img.shields.io/badge/preprint-arXiv-gold.svg)](https://arxiv.org/abs/2105.13469)

---

## Overview

DTAmc is an R package to perform statistical analysis for Diagnostic Test Accuracy
(DTA) studies with co-primary endpoints sensitivity and specificity. 
In particular, it allows for multiplicity corrections when the diagnostic accuracy of multiple index tests (or candidate prediction models) is investigated simultaneously. 
The implemented methods are also applicable in the extended case of more than
two subpopulations of interest (i.e. multiclass classification problems).

**NOTE: THIS IS AN EARLY DEVELOPMENT VERSION. VALIDITY OF RESULTS IS NOT GUARANTEED!**


---

## Installation

This package can be installed from GitHub with the following command

```r
remotes::install_github('maxwestphal/DTAmc', build_vignettes = TRUE)
```

---

## Usage

This package contains a vignette which explains the basic functionality.
The vignette can be displayed with the command

```r
vignette("DTAmc")
```

---

## References

- [Westphal, Max, and Antonia Zapf. "Statistical Inference for Diagnostic Test Accuracy Studies with Multiple Comparisons." arXiv preprint arXiv:2105.13469 (2021).](https://arxiv.org/abs/2105.13469)
