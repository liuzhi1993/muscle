# muscle verification

## Correctness

The package is compared against the original `muscle` source package rather
than against duplicated expected values.

- All exported helper, simulation, `MUSCLE`, split, full-interval,
  deconvolution, and `MMUSCLE` tests pass.
- PST range-quantile queries are exhaustively compared with brute force for
  prefix, shifted-one, and shifted-two windows.
- The exhaustive tests cover every valid order, every query endpoint pair,
  both start parities, repeated observations, maximum queries, minimum
  queries, and leftmost tie-breaking.
- Additional reference comparisons cover constant, tied, piecewise, outlier,
  and randomized signals.
- `R CMD check --no-manual muscle_0.1.1.tar.gz` completes with `Status: OK`.

The original shifted-two wavelet-tree path has undefined behavior when an
extreme quantile requests order `l` from a window containing only `l - 1`
observations. It can terminate R with a segmentation fault. `muscle`
defines this boundary deterministically as the window maximum. Strict output
comparison is therefore performed where the original routine has a defined
result.

## Gaussian benchmark

Configuration:

- `n = 2000`
- change points at 986 and 1016
- levels `-4`, `0`, and `4`
- Gaussian noise variance `0.9`
- `alpha = 0.3`, `beta = 0.5`
- seed `20260726`

Observed elapsed times from one reproducible run:

| Backend | Time (seconds) | Speed relative to original |
|---|---:|---:|
| `muscle` (`original`) | 14.307 | 1.0x |
| `muscle` (`PRT`) | 0.101 | 141.7x |
| `muscle` (`PST`) | 0.731 | 19.6x |

For this run, the `PST` backend and the original implementation have identical output names
and `left`, and equal `value` within tolerance `1e-12`.

At `n = 2000`, the existing nonpersistent segment-tree cache is faster. The
PST implementation is intended to reduce worst-case cache space from
`O(n^2)` to `O(n log^2 n)`, trading an additional logarithmic query factor for
substantially better asymptotic storage.
