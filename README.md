# IPFStructPenalty

`IPFStructPenalty` fits the Lasso, elastic and tree-lasso methods with penalty factors for different data sources, and uses the Efficient Parameter Selection via Global Optimization (EPSGO) algorithm to optimize multiple penalty parameters.

## R package instruction

You can download the package and install it locally, or install it from github directly by the code:

```{r setup, include=FALSE}
library(devtools)
devtools::install_github("zhizuio/IPFStructPenalty/IPFStructPenalty")
```
See the file `test.R` for s simple start with simulated data. The scripts `Sim_reg.R` and `GDSC_reg.R` are for simulation study and analyzing real data, respectively, presented in [Zhao \& Zucknick (2019)](https://arxiv.org/abs/1902.04996)..

## Reference
The pre-print has been arXived [Zhao \& Zucknick (2019)](https://arxiv.org/abs/1902.04996): https://arxiv.org/abs/1902.04996.

## Update
### New in the version IPFStructPenalty_0.1.2.tar.gz (02 August 2019):

The help files of functions have been updated.

### New in the version IPFStructPenalty_0.1.1.tar.gz (01 August 2019):

1) The `DESCRIPTION` file has been updated.

2) The conditional logistic lasso has been added, but using the equivalent Cox lasso model because the stratification factor.

