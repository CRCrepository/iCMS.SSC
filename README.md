# iCMS single-sample classifier (SSC)

The single-sample classifier (SSC) estimates the most probable intrinsic consensus molecular subtype (iCMS) of a colorectal cancer sample as described in Joanito et al. (Nat Genet. 2022 Jul;54(7):963-975). 

The package has been developed to be relatively robust to technical variation, specifically without reliance on batch correction methods. An implementation of the bulk iCMS classifier is also provided, based on the NTP method, as described in the original publication.


## installation

Please follow the instruction below to install iCMS single-sample classifier.

```{r}
if(!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github('CRCrepository/iCMS.SSC')
```

## Tutorial

The algorithm for the estimation of the iCMS is based on distance from a reference set synthetic exemplar centroids. There are two criteria that can be applied, based on correlation distance quantiles (DQ) or K-nearest neighbors voting (KNN). The KNN criterion is slightly faster than DQ. The published bulk NTP algorithm is also provided.

All three functions accept gene expression matrices with log-transformed measures of expression. Genes should be in rows (gene symbol row names) and samples in columns. Note that the published bulk NTP algorithm is relatively sensitive to batch effect. Appropriate measures should be considered if you have reason to suspect that the input matrix is influenced by batch effect.
Please see the documentation of each individual function for additional details. 

Example code below:

```{r}
library(iCMS.SSC)

data(test) ## load test dataset
icms_results.dq = iCMS.DQ(test)
icms_results.knn = iCMS.KNN(test)
icms_results.ntp = iCMS.NTP(test)

```
The outputs provide iCMS calls for each algorithm (DQ, KNN, and
NTP). All three functions provide both nearest (most likely) and
confident iCMS calls in the output under the columns "nearest.icms"
and "confident.icms" respectively. The difference between the nearest
and confident calls is that a nearest call will always be made, even
when iCMS2 and iCMS3 classes seem almost equally probable. A confident
call is only returned if it seems statistically quite probable. When a
confident call cannot be made, the "confident.icms" column will
contain NA (not a value). The "nearest.icms" column should not contain
NAs.

A sensible choice in most cases is to call iCMS.KNN() with default
parameters. The iCMS.DQ() function produces fewer confident calls than
the iCMS.KNN (~80% versus 90+%), but has a slightly higher accuracy
when a confident call is made (~99% vs 95%). 

Note that the output data frames also contain variables that are
specific to each algorithm/criterion, mostly useful for debugging and
deeper analysis. These can be safely ignored in most cases.


```{r}
icms_results.dq$nearest.icms ## the most probable iCMS class, including low confidence
icms_results.dq$confident.icms ## high-confidence iCMS class
``` 

For more details please see the individual documentation of each function in R. 
