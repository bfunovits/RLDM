# Blanchard/Quah (1989) dataset

This dataset contains a (159 x 2)-dimensional matrix with quarterly data
of real GDP growth rates (first column) and the unemployment rate in the
USA. It starts in the second quarter of 1948 and ends in the last
quarter of 1987. The script for transforming the raw data (as csv-file)
to the matrix BQdata is available in the *data-raw* directory. This
dataset was used in [Gourieroux, Monfort, Renne
(2019)](https://doi.org/10.1093/restud/rdz028)

## Usage

``` r
BQdata

BQdata_ts

BQdata_xts
```

## Format

A tibble (BQdata), a ts-object (BQdata_ts), or an xts-object
(BQdata_xts) with 159 rows and 2 variables (plus a timestamp, so in
total 3):

- date:

  Timestamp as column of type *dttm* for the tibble, and as index for
  the *ts* and *xts* object.

- rGDPgrowth_demeaned:

  The real GDP growth series is demeaned with respect to two different
  subperiods: Till the last quarter of 1973 and from the first quarter
  of 1974

- unemp_detrended:

  The unemployment rate is detrended with respect to a linear trend (and
  an intercept).

An object of class `mts` (inherits from `ts`, `matrix`) with 159 rows
and 2 columns.

An object of class `xts` (inherits from `zoo`) with 159 rows and 2
columns.

## Source

<https://academic.oup.com/restud/advance-article/doi/10.1093/restud/rdz028/5490841#supplementary-data>
