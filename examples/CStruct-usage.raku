#!/usr/bin/env raku
use v6.d;

use lib <. lib>;
use NativeCall;
use Math::SparseMatrix::Native;

my $matrix1 = Math::SparseMatrix::Native::CSRStruct.new().random(nrow => 1000, ncol => 10_000, nnz => 20000);

say (nrow => $matrix1.nrow, ncol => $matrix1.ncol, nnz => $matrix1.nnz);
#$matrix4.print;

my $tstart = now;
my $matrix1t = $matrix1.transpose;
my $tend = now;
say (:$matrix1t);
say (nrow => $matrix1t.nrow, ncol => $matrix1t.ncol, nnz => $matrix1t.nnz);
say "Transpose time: {$tend - $tstart} seconds.";


$tstart = now;
my $matrix2 = $matrix1.dot-pattern($matrix1t);
$tend = now;
say (:$matrix2);
say (nrow => $matrix2.nrow, ncol => $matrix2.ncol, nnz => $matrix2.nnz);
say "Dot-pattern time: {$tend - $tstart} seconds.";

$tstart = now;
my $matrix3 = $matrix1.dot-pattern($matrix1t);
$tend = now;
say (:$matrix3);
say (nrow => $matrix3.nrow, ncol => $matrix3.ncol, nnz => $matrix3.nnz);
say "Dot-numeric time: {$tend - $tstart} seconds.";