#!/usr/bin/env raku
use v6.d;

# Using tolerance => 10e-12 : ≈ 6 times slower that the corresponding Wolfram Language profiling results.
# Using tolerance => 10e-8 : ≈ 6 times slower that the corresponding Wolfram Language profiling results.
#`[

{m, n} = {1000, 600};
density = 0.2;
nnz = Round[m*n*density];
triplets = Transpose@{RandomInteger[{1, n}, nnz], RandomInteger[{1, m}, nnz], RandomReal[{0, 1}, nnz]};
matrix = SparseArray[Map[Most[#] -> Last[#] &, triplets]];
MatrixForm[matrix[[1 ;; 20, 1 ;; 20]]]

Length@matrix["NonzeroValues"]

(* Out[136]= 108759 *)

AbsoluteTiming[
  res = SingularValueDecomposition[matrix, 100];
]
(* {2.83412, Null} first time call*)
(* {1.99494, Null} subsequent times *)

]

#use lib <. lib>;
use Math::SparseMatrix::Native;

my $nrow = 1000;
my $ncol = 600;
my $density = 0.2;
my $nnz = ($nrow * $ncol * $density).Int;
my $seed = 3432;
my $k = 100;
my $tolerance = 1e-16;

my $tstart = now;
my $matrix1 = Math::SparseMatrix::Native::CSRStruct.new.random(:$nrow, :$ncol, :$nnz, :$seed);
my $tend = now;
say (:$matrix1);
say "Creation time: { $tend - $tstart } seconds.";
say "Non-zero values 1: ", $matrix1.explicit-length;
say "Fill in 1: ", $matrix1.explicit-length / $matrix1.rows-count / $matrix1.columns-count;
say "-" x 100;

$tstart = now;
my ($u, $s, $v) = $matrix1.svd($k, :$tolerance);
$tend = now;
say "SVD time: { $tend - $tstart } seconds.";


.say for |(:$u, :$s, :$v);
