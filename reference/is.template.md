# Check templates

Check templates

## Usage

``` r
is.template(tmpl, class = c("any", "stspmod", "armamod", "rmfdmod"))
```

## Arguments

- tmpl:

  object to be tested

- class:

  test for a specific class

## Value

TRUE/FALSE

## Examples

``` r
is.template(1)
#> [1] FALSE
is.template(tmpl_llm(), 'armamod')
#> [1] FALSE
is.template(tmpl_llm(), 'stspmod')
#> [1] TRUE
```
