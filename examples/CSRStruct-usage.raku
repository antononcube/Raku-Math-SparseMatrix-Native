#!/usr/bin/env raku
use v6.d;

use lib <. lib>;
use NativeCall;
use Math::SparseMatrix::CSR;
use Math::SparseMatrix::Native;
use NativeHelpers::Array;

sub normalize(Math::SparseMatrix::Native::CSRStruct:D $matrix1) {
   return Math::SparseMatrix::CSR.new(
            row-ptr => copy-to-array($matrix1.row_ptr, $matrix1.nrow + 1)>>.Num,
            col-index => copy-to-array($matrix1.col_index, $matrix1.nnz),
            values => copy-to-array($matrix1.values, $matrix1.nnz),
            nrow => $matrix1.nrow,
            ncol => $matrix1.ncol,
            implicit-value => $matrix1.implicit_value
            );
}


my $matrix1 = Math::SparseMatrix::Native::CSRStruct.new().random(nrow => 1000, ncol => 10_000, nnz => 20_000);
#my $matrix1 = Math::SparseMatrix::Native::CSRStruct.new().random(nrow => 6, ncol => 10, nnz => 10);
say (nrow => $matrix1.nrow, ncol => $matrix1.ncol, nnz => $matrix1.nnz);

if $matrix1.nnz < 100 { $matrix1.&normalize.print }

say '-' x 100;

my $tstart = now;
my $matrix1t = $matrix1.transpose;
my $tend = now;
say (:$matrix1t);
say (nrow => $matrix1t.nrow, ncol => $matrix1t.ncol, nnz => $matrix1t.nnz);
say "Transpose time: {$tend - $tstart} seconds.";

if $matrix1t.nnz < 100  { $matrix1t.&normalize.print }

say '-' x 100;

$tstart = now;
my $matrix2 = $matrix1.dot-pattern($matrix1t);
$tend = now;
say (:$matrix2);
say (nrow => $matrix2.nrow, ncol => $matrix2.ncol, nnz => $matrix2.nnz);
say "Dot-pattern time: {$tend - $tstart} seconds.";

if $matrix2.nnz < 100 {
    $matrix2.&normalize.print
}

say '-' x 100;

$tstart = now;
my $matrix3 = $matrix1.dot($matrix1t);
$tend = now;
say (:$matrix3);
say (nrow => $matrix3.nrow, ncol => $matrix3.ncol, nnz => $matrix3.nnz);
say "Dot-numeric time: {$tend - $tstart} seconds.";

if $matrix3.nnz < 100 {
    $matrix3.&normalize.print
}

say '-' x 100;

my $matrix4 = $matrix1.add(100, :clone);
say (:$matrix4);
say (implicit_value => $matrix4.implicit_value);

if $matrix4.nnz < 100 {
    $matrix4.&normalize.print(:!iv)
}

say '-' x 100;

my $matrix5 = $matrix1.add($matrix1);
say (:$matrix5);

if $matrix5.nnz < 100 {
    $matrix5.&normalize.print(:!iv)
}