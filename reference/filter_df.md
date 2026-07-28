# Filter a long-format data frame to valid spatial cells

Keeps only rows whose spatial index (`ids`) appears in `valid_ids` and
adds a `sid` column (1..n_valid) giving each valid cell's rank among the
valid cells. `sid` matches the indexing of `membership` vectors returned
by [`genclust()`](genclust.md). If `valid_ids` is `NULL`, all rows are
kept and `sid` is assigned by rank of `ids`.

## Usage

``` r
filter_df(df, valid_ids = NULL)
```

## Arguments

- df:

  A long-format data frame as returned by [`data_all()`](data_all.md).

- valid_ids:

  Integer vector of valid flat spatial positions, as returned in
  `genclust(...)$valid_ids`. If `NULL`, all spatial cells are treated as
  valid.

## Value

A filtered data frame with rows corresponding to valid spatial cells
only, with an additional `sid` column (integer, 1..n_valid).
