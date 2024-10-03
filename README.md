# Math::SparseMatrix::Native

[![Actions Status](https://github.com/antononcube/Raku-SparseMatrix-Native/actions/workflows/linux.yml/badge.svg)](https://github.com/antononcube/Raku-SparseMatrix-Native/actions)
[![Actions Status](https://github.com/antononcube/Raku-SparseMatrix-Native/actions/workflows/macos.yml/badge.svg)](https://github.com/antononcube/Raku-SparseMatrix-Native/actions)

[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic%202.0-0298c3.svg)](https://opensource.org/licenses/Artistic-2.0)

Raku package with sparse matrix algebra functions implemented in C.
Apple's Accelerate library is used if available.

------

## TODO

- [ ] TODO Core functionalities
    - [X] DONE C-struct representation: data and methods
    - [X] DONE `transpose`
    - [X] DONE `dot-pattern`
    - [ ] TODO `dot` matrix x dense vector
    - [X] DONE `dot` (and `dot-numeric`) matrix x matrix
    - [X] DONE `add` with a scalar
    - [X] DONE `add` with another matrix
    - [ ] TODO `multiply` with a scalar
    - [ ] TODO `multiply` with another matrix
- [ ] TODO Adaptation to "Math::SparseMatrix"
    - This package was made in order to have faster computation with "Math::SparseMatrix".
    - But it can be self-contained and independent from "Math::SparseMatrix".
    - Hence, we make an adapter class in "Math::SparseMatrix".
- [ ] TODO Unit tests
    - [ ] TODO Creation (and destruction)
    - [ ] TODO Element-wise operations
    - [ ] TODO Dot product

------

## References

[AAp1] Anton Antonov,
[Math::SparseMatrix Raku package](https://github.com/antononcube/Raku-Math-SparseMatrix),
(2024),
[GitHub/antononcube](https://github.com/antononcube).





