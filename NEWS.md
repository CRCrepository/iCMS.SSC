# iCMS.SSC 1.2

## New features

* `iCMS.DQ()` / `dq.vect()`: new `stratify.ds` parameter — computes the
  quantile margin separately for each training dataset (Marisa, PETACC, TCGA)
  and aggregates, removing confounding by platform-level similarity differences.
* `iCMS.DQ()` / `dq.vect()`: new `weighted.genes` parameter — uses
  weighted Pearson or weighted cosine similarity (all genes, weighted by
  absolute t-statistic) for `metric = "pearson"` / `"cosine"`; restricts
  to the top 100 most discriminative genes for rank-based metrics.
* `iCMS.DQ()` / `dq.vect()`: `verbose` parameter now defaults to `FALSE`,
  returning 7 summary columns instead of the full 192-column prototype matrix.
  Set `verbose = TRUE` to restore the previous behaviour.
* `iCMS.KNN()` / `knn.vect()`: new `stratify.ds` parameter — each training
  dataset casts one vote (best i2 vs. best i3 prototype); unanimous agreement
  across all three datasets produces a confident call.
* `iCMS.KNN()` / `knn.vect()`: new `weighted.genes` parameter (same
  semantics as the DQ method above).
* Added internal `gene.weights` data object: named vector of absolute Welch
  t-statistics for all 201 discriminative genes, used by `weighted.genes`.
* Added internal helper functions `.cosine_vec_mat()`, `.wt_pearson_mat()`,
  `.wt_cosine_mat()` for efficient vectorised similarity computation.

## Improvements

* `iCMS.KNN()` / `knn.vect()`: neighbours are now **distance-weighted** —
  each of the `nn` neighbours contributes a vote proportional to its
  (non-negative shifted) similarity score instead of equal weight.
* `iCMS.KNN()` / `knn.vect()`: confident call criterion changed from a
  hardcoded 90% majority threshold to a **binomial test**
  (`pbinom(n_majority - 1, nn, 0.5, lower.tail = FALSE) < 0.05`). This is
  statistically motivated and consistent across different values of `nn`.

## Bug fixes

* `iCMS.KNN()`: the `verbose` argument was silently routed to the `mcenter`
  positional argument in the internal dispatch call, causing `verbose = TRUE`
  to inadvertently mean-centre samples. Fixed by switching to fully named
  argument passing.

# iCMS.SSC 1.1

* Fixed cosine distance calculation in `calc_knn.R`.
* Corrected cosine similarity for vector inputs.
