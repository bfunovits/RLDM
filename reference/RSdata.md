# Ramey/Shapiro dataset

This dataset contains a (248 x 7)-dimensional matrix with quarterly US
data from the first quarter of 1947 till the last quarter of 2008. The
script for transforming the raw data (as csv-file) to the matrix RSdata
is available in the *data-raw* directory. The data set is used to
identify government spending shocks and in particular to investigate the
response of consumption and real wages with respect to these shocks. All
variables are orthogonalized with respect to an intercept and a linear
trend. See <https://doi.org/10.1093/qje/qjq008> for more details.

## Usage

``` r
RSdata

RSdata_ts

RSdata_xts
```

## Format

A matrix with 248 rows and 7 variables (plus a timestamp, so in total 8
variables). Versions as tibble (RSdata), ts (RSdata_ts), or xts
(RSdata_xts) object are available:

1.  Date column

2.  the logarithm of real per capita quantities of total government
    spending

3.  real GDP

4.  total hours worked

5.  nondurable plus services consumption

6.  private fixed investment

7.  real wage (or more precisely the nominal compensation in private
    business divided by the deflator in private business)

8.  tax rate

An object of class `mts` (inherits from `ts`, `matrix`) with 248 rows
and 7 columns.

An object of class `xts` (inherits from `zoo`) with 248 rows and 7
columns.

## Source

<https://econweb.ucsd.edu/~vramey/research.html>
