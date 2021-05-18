## DTAmc 

DTAmc is an R package to perform statistical analysis for Diagnostic Test Accuracy
(DTA) studies with co-primary endpoints sensitivity and specificity. 
In particular,
it allows for multiplicity corrections when the diagnostic accuracy of multiple 
index tests (or candidate models) is investigated simultaneously. 
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

- Max Westphal, Antonia Zapf, 2021, Statistical Inference for Diagnostic Test Accuracy Studies with Multiple Comparisons, in preparation.