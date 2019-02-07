# IPFStructPenalty

`IPFStructPenalty` fits Lasso, elastic and tree-lasso methods with penalty factors for different data sources, and uses the Efficient Parameter Selection via Global Optimization (EPSGO) algorithm to optimize multiple penalty parameters.

## R package instruction

You can download and install this package locally, or install it from github directly with:

```{r setup, include=FALSE}
devtools::install_github("zhizuio/IPFStructPenalty/IPFStructPenalty")
```
See the file `Test.R` for usage with simulated data, and `GDSC_reg.R` for analyzing real data.
