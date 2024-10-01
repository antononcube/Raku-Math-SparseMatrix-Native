#!/usr/bin/env raku
use v6.d;

#use lib <. lib>;
use Math::SparseMatrix::Native;
use Math::SparseMatrix::Utilities;
use NativeCall;
use NativeHelpers::Array;
use JSON::Fast;

my $nrow = 1_000;
my $ncol = 10_000;
my $density = 0.001;
my $tol = 0.01;
my $type = 'CSR';

say "-" x 100;
say "Matrix 1:";
my $tstart = now;
my $matrix1 = generate-random-sparse-matrix($nrow, $ncol, :$density, :$tol, :$type):!decorated;
say $matrix1;
my $tend = now;
say "Generation time: {$tend - $tstart} seconds.";
say "-" x 100;
$matrix1.print if $matrix1.ncol < 20 && $matrix1.nrow < 20;
say "-" x 100;

$tstart = now;
my $matrix2 = Math::SparseMatrix::Native::CSR.new(
        row-ptr => $matrix1.row-ptr,
        col-index => $matrix1.col-index,
        values => $matrix1.values,
        nrow => $matrix1.nrow,
        ncol => $matrix1.ncol,
        implicit-value => $matrix1.implicit-value
        );
$tend = now;
say $matrix2.raku;
say "Creation time: {$tend - $tstart} seconds.";

$tstart = now;
my $matrix3 = $matrix2.transpose;
$tend = now;
say (:$matrix3);
say "Transpose time: {$tend - $tstart} seconds.";


say "-" x 100;

my $matrix4 = Math::SparseMatrix::CSR.new(
        row-ptr => $matrix3.row-ptr.Array,
        col-index => $matrix3.col-index.Array,
        values => $matrix3.values.Array,
        nrow => $matrix3.nrow,
        ncol => $matrix3.ncol,
        implicit-value => $matrix3.implicit-value
        );

say $matrix4;
$matrix4.print if $matrix4.ncol < 20 && $matrix4.nrow < 20;

say "-" x 100;

$tstart = now;
my $matrix5 = $matrix2.dot-pattern($matrix3);
$tend = now;
say (:$matrix5);
say "Dot-pattern time: {$tend - $tstart} seconds.";

my $matrix6 = Math::SparseMatrix::CSR.new(
        row-ptr => $matrix5.row-ptr.Array,
        col-index => $matrix5.col-index.Array,
        values => $matrix5.values.Array,
        nrow => $matrix5.nrow,
        ncol => $matrix5.ncol,
        implicit-value => $matrix5.implicit-value
        );

$matrix6.print if $matrix5.ncol < 20 && $matrix5.nrow < 20;

say "-" x 100;
