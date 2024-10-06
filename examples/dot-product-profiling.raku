use Math::SparseMatrix::Native;
use Math::SparseMatrix::Native::Utilities;

my $nrow = 1000;
my $ncol = 10_000;
my $density = 0.005;
my $nnz = ($nrow * $ncol * $density).Int;
my $seed = 3432;
my $n = 1000;

my $matrix1 = Math::SparseMatrix::Native::CSRStruct.new.random(:$nrow, :$ncol, :$nnz, :$seed);
my $matrix2 = Math::SparseMatrix::Native::CSRStruct.new.random(nrow => $ncol, ncol => $nrow, :$nnz, :$seed);

say (:$matrix1);
say (:$matrix2);


my $tstart = now;
for ^$n {
    $matrix1.dot($matrix2)
}
my $tend = now;

say "Total time : { $tend - $tstart }";
say "Mean time  : { ($tend - $tstart) / $n }"