# Math::SparseMatrix::Native

[![Actions Status](https://github.com/antononcube/Raku-Math-SparseMatrix-Native/actions/workflows/linux.yml/badge.svg)](https://github.com/antononcube/Raku-Math-SparseMatrix-Native/actions)
[![Actions Status](https://github.com/antononcube/Raku-Math-SparseMatrix-Native/actions/workflows/macos.yml/badge.svg)](https://github.com/antononcube/Raku-Math-SparseMatrix-Native/actions)

[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic%202.0-0298c3.svg)](https://opensource.org/licenses/Artistic-2.0)

Raku package with sparse matrix algebra functions implemented in C.

The package is a "standalone" implementation of sparse matrix functionalities using 
Compressed Sparse Row (CSR) format. The intent, though, is to use the class `Math::SparseMatrix::Native::CSRStruct`
provided by this package in the general package 
["Math::SparseMatrix"](https://github.com/antononcube/Raku-Math-SparseMatrix), [AAp1].
(See the section "Design" below.)

The core algorithms of the package follow the FORTRAN implementations given in the book
"Sparse Matrix Technology" by Pissanetzky, [SP1].

This package should be compared -- and likely replaced -- with Raku package(s) interfacing 
sparse matrix algebra libraries, like, 
["SuiteSparse"](https://github.com/DrTimothyAldenDavis/SuiteSparse).
(When/if such Raku packages are implemented.)

**Remark:** This package uses a C implementation based on standard C libraries ("math", "stdio", "stdlib", "string".)
Hence, it should be cross platform.

**Remark:** Currently, on macOS, Apple's Accelerate library is _not_ used.
There are several reasons for this: 
(i) lack of appropriate documentation to sparse linear algebra in C,
(i) using dense matrices for sparse matrix computations with its documented, older LAPACK interface libraries.

------

## Usage examples

Make two random -- relatively large -- sparse matrices:

```perl6
use Math::SparseMatrix::Native;
use Math::SparseMatrix::Native::Utilities;

my $nrow = 1000;
my $ncol = 10_000;
my $density = 0.005;
my $nnz = ($nrow * $ncol * $density).Int;
my $seed = 3432;

my $matrix1 = Math::SparseMatrix::Native::CSRStruct.new.random(:$nrow, :$ncol, :$nnz, :$seed);
my $matrix2 = Math::SparseMatrix::Native::CSRStruct.new.random(nrow => $ncol, ncol => $nrow, :$nnz, :$seed);

say (:$matrix1);
say (:$matrix2);
```
```
# matrix1 => Math::SparseMatrix::Native::CSRStruct(:specified-elements(50000), :dimensions((1000, 10000)), :density(0.005))
# matrix2 => Math::SparseMatrix::Native::CSRStruct(:specified-elements(50000), :dimensions((10000, 1000)), :density(0.005))
```

**Remark:** Compare the dimensions, densities, and the number of "specified elements" in the gists. 

Here are 100 dot-products of those matrices (with timings):

```perl6
my $tstart=now;
my $n = 100;
for ^$n {
    $matrix1.dot($matrix2)
}
my $tend=now;
say "Total time : {$tend - $tstart}";
say "Mean time  : {($tend - $tstart)/$n}"
```
```
# Total time : 0.187534002
# Mean time  : 0.00187534002
```

------

## Design

The primary motivation for implementing this package is the need to have _fast_ sparse matrix computations.

- The pure-Raku implemented sparse matrix algorithms in 
["Math::SparseMatrix"](https://raku.land/zef:antononcube/Math::SparseMatrix), [AAp1],
are too slow. 

- Same algorithms implemented in C (in this package) are ≈100 times faster.  

The `NativeCall` class `Math::SparseMatrix::Native::CSRStruct` has Compressed Sparse Row (CSR) format sparse matrix operations.
The [Adapter pattern](https://en.wikipedia.org/wiki/Adapter_pattern) is used to include `Math::SparseMatrix::Native::CSRStruct`
into the `Math::SparseMatrix` hierarchy:

```mermaid
classDiagram
    class Abstract["Math::SparseMatrix::Abstract"] {
        <<abstract>>
        +value-at()
        +row-at()
        +column-at()
        +row-slice()
        +AT-POS()
        +print()
        +transpose()
        +add()
        +multiply()
        +dot()
    }

    class CSRStruct {
        <<C struct>>
    }
    
    class NativeCSR["Math::SparseMatrix::CSR::Native"] {
        $row_ptr
        $col_index
        @values
        nrow
        ncol
        implicit_value
    }

    class NativeAdapater["Math::SparseMatrix::NativeAdapter"] {
        +row-ptr
        +col-index
        +values
        +nrow
        +ncol
        +implicit-value
    }    
    
    class SparseMatrix["Math::SparseMatrix"] {
        Abstract core-matrix
        +AT-POS()
        +print()
        +transpose()
        +add()
        +multiply()
        +dot()
    }
    
    NativeAdapater --> Abstract : implements
    SparseMatrix --> Abstract : Hides actual component class
    SparseMatrix *--> Abstract
    NativeAdapater *--> NativeCSR
    NativeCSR -- CSRStruct : reprents
```

------

## TODO

- [X] DONE Core functionalities
    - [X] DONE C-struct representation: data and methods
    - [X] DONE `transpose`
    - [X] DONE `dot-pattern`
    - [X] DONE `dot` matrix x dense vector
    - [X] DONE `dot` (and `dot-numeric`) matrix x matrix
    - [X] DONE `add` with a scalar
    - [X] DONE `add` with another matrix
    - [X] DONE `multiply` with a scalar
    - [X] DONE `multiply` with another matrix
    - [X] DONE Info methods
    - [X] DONE Access functions
    - [X] DONE Values operations `unitize`, `clip`, `round`.
- [X] DONE Refactoring
  - Consistent use of `unsigned int` or `int` for `row_ptr` and `col_index`. 
- [X] DONE Adaptation to "Math::SparseMatrix"
    - This package was made in order to have faster computation with "Math::SparseMatrix", [AAp1].
    - But it can be self-contained and independent from "Math::SparseMatrix".
    - Hence, we make an adapter class:
      - See `Math::SparseMatrix::NativeAdapter` in "Math::SparseMatrix".
- [X] DONE Unit tests
    - [X] DONE Creation (and destruction)
    - [X] DONE Access
    - [X] DONE Element-wise operations
    - [X] DONE Dot product
    - [X] DONE Values operations `unitize`, `clip`, `round`.

------

## References

### Books

[SP1] Sergio Pissanetzky, Sparse Matrix Technology, Academic Pr (January 1, 1984), ISBN-10: 0125575807, ISBN-13: 978-0125575805.

### Packages

[AAp1] Anton Antonov,
[Math::SparseMatrix Raku package](https://github.com/antononcube/Raku-Math-SparseMatrix),
(2024),
[GitHub/antononcube](https://github.com/antononcube).

