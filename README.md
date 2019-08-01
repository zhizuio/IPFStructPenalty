# IPFStructPenalty

`IPFStructPenalty` fits the Lasso, elastic and tree-lasso methods with penalty factors for different data sources, and uses the Efficient Parameter Selection via Global Optimization (EPSGO) algorithm to optimize multiple penalty parameters.

## R package instruction

You can download the package and install it locally, or install it from github directly by the code:

```{r setup, include=FALSE}
library(devtools)
devtools::install_github("zhizuio/IPFStructPenalty/IPFStructPenalty")
```
See the file `Example.R` for s simple start with simulated data. The scripts `Sim_reg.R` and `GDSC_reg.R` are for simulation study and analyzing real data, respectively, presented in [Zhao \& Zucknick (2019)](https://arxiv.org/abs/1902.04996)..

## Reference
The pre-print has been arXived [Zhao \& Zucknick (2019)](https://arxiv.org/abs/1902.04996).
