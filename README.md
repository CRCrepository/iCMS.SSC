# iCMS single-sample classifier (SSC)

Need to explain

## installation

Please follow the instruction below to install iCMS single-sample classifier.

```{r}
if(!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github('CRCrepository/iCMS.SSC')
```

## Tutorial

Follow these instructions to use the examples in the package.
The input should be a normalised and log-transformed matrix, having with gene symbol (row) and sample name (column).

```{r}
library(iCMS.SSC)

data(test) ## load test dataset
icms_results.dq = iCMS.DQ(test)
icms_results.knn =iCMS.KNN(test)
icms_results.ntp = iCMS.NTP(test)

```
The outputs provide iCMS calls for each algorithm (DQ, KNN, and NTP). 
