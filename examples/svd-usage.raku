#!/usr/bin/env raku
use v6.d;

use Math::SparseMatrix::Native;

my @m = [0.5,0.5,0.5,0,0],[1.0,2.0,4.0,8.0,0],[1.5,4.5,13.5,40.5,121.5],[0,8.0,32.0,128.0,512.0],[0,0,62.5,312.5,1562.5];

my $matrix = Math::SparseMatrix::Native::CSRStruct.new(dense-matrix => @m);
my ($u, $s, $v) = $matrix.svd(3);
my $k = min($matrix.nrow, $matrix.ncol);


my $tol = 0.0001;

say (:$u);
.say for |$u.Array».round($tol);

say (:$s);
.say for |$s.Array».round($tol);

say (:$u);
.say for |$u.Array».round($tol);
