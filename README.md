# IPFStructPenalty

`IPFStructPenalty` fits the Lasso, elastic and tree-lasso methods with penalty factors for different data sources, and uses the Efficient Parameter Selection via Global Optimization (EPSGO) algorithm to optimize multiple penalty parameters.

## R package instruction

You can download the package and install it locally, or install it from github directly by the code:

```{r setup, include=FALSE}
library(devtools)
devtools::install_github("zhizuio/IPFStructPenalty/IPFStructPenalty")
```
See the file `test.R` for s simple start with simulated data. The scripts `Sim_reg.R` and `GDSC_reg.R` are for simulation study and analyzing real data, respectively, presented in [Zhao \& Zucknick (2020)](https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssc.12400). The original runnings were on HPC. They will be every slow to be run on PC.

## Reference
Zhao Z, Zucknick M (2020). *Structured penalized regression for drug sensitivity prediction*. Journal of the Royal Statistical Society: Series C (Applied Statistics). 69(3): 525-545, URL: https://doi.org/10.1111/rssc.12400.

## Update
### New in version IPFStructPenalty_1.0.2.tar.gz (15 November 2019):

Fixed the specified directory of the installed package (`lib.loc`) for the parallelization.

### version IPFStructPenalty_1.0.1.tar.gz (30 August 2019):

The unpenalized argument for the conditional logistic lasso has been added and passed from the EPSGO algorithm.

### Version IPFStructPenalty_1.0.tar.gz (23 August 2019):

The package and all codes have been adapted to do all data analysis in [Zhao \& Zucknick (2020)](https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssc.12400).

### Version IPFStructPenalty_0.1.4.tar.gz (21 August 2019):

The parallelisation issue of the conditional logistic lasso have been fixed.

### Version IPFStructPenalty_0.1.3.tar.gz (15 August 2019):

Some bugs of the conditional logistic lasso have been fixed.

### Version IPFStructPenalty_0.1.2.tar.gz (02 August 2019):

The help files of functions have been updated.

### Version IPFStructPenalty_0.1.1.tar.gz (01 August 2019):

1) The `DESCRIPTION` file has been updated.

2) The conditional logistic lasso has been added, but using the equivalent Cox lasso model because the stratification factor.

