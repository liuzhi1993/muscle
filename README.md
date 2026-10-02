# muscle 0.1.3

This package exposes one `MUSCLE` interface with three selectable backends.
The default is `implementation = "PRT"`, which uses the parity range tree
cache. `implementation = "PST"` selects the persistent segment tree backend,
and `implementation = "original"` selects the original wavelet-tree search.

The existing `dyadic`, `split`, `deconv`, and `MMUSCLE` interfaces remain
available for every backend.

For each dyadic window length, versions are created by activating observations
in descending value order. Activating observation `Y[t]` performs one range-add
on the contiguous set of window starts whose windows contain `t`. Unchanged
subtrees are shared between versions. Each node stores lazy range-add state and
the minimum and maximum active counts for even and odd window starts.

A range quantile query binary-searches the threshold versions and performs a
range extrema query at each step. Ties are grouped into one version, and the
leftmost window start is returned when several windows have the same extremal
quantile, matching the original implementation.

For dyadic lengths, the PST cache has worst-case construction time and space
`O(n log^2 n)`. A range quantile query takes `O(log^2 n)`, and one dyadic
multiscale test on a fixed candidate interval takes `O(log^3 n)`.

The test suite includes exhaustive PST-versus-brute-force checks and reference
comparisons against the original package for every exported function. The
original shifted-two implementation has undefined behavior when it requests
order `l` from a window of length `l - 1` at extreme quantile levels;
`muscle` handles that boundary deterministically as the window maximum.
