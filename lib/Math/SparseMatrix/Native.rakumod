use v6.d;

unit module Math::SparseMatrix::Native;

use NativeCall;
use NativeHelpers::Array;

my constant $library = %?RESOURCES<libraries/SparseMatrixFunctions>;

# void dot_pattern(int *IA, int *JA, int *IB, int *JB, int NP, int NQ, int NR, int *IC, int *JC)
sub dot-pattern(CArray[int32], CArray[int32], CArray[int32], CArray[int32],
                int32, int32, int32,
                CArray[int32], CArray[num64]) is native($library) is symbol('dot_pattern') {*}

#sub DotProduct(CArray[num64], CArray[num64], int32 --> num64) is native($library) {*}
#sub DotProductFloat(CArray[num32], CArray[num32], int32 --> num32) is native($library) {*}

sub transpose(CArray[int32], CArray[int32], CArray[num64],
                uint32, uint32,
                CArray[int32], CArray[int32], CArray[num64]) is native($library) {*};


#-----------------------------------------------------------

class CSR {
    has CArray[int32] $.row-ptr;
    has CArray[int32] $.col-index;
    has CArray[num64] $.values;
    has UInt:D $.nrow is required;
    has UInt:D $.ncol is required;
    has Numeric:D $.implicit-value = 0.0;

    method new(:@values, :@col-index, :@row-ptr, UInt:D :$nrow, UInt:D :$ncol, Numeric:D :$implicit-value = 0.0) {
        self.bless(
                row-ptr => CArray[int32].new(@row-ptr),
                col-index => CArray[int32].new(@col-index),
                values => CArray[num64].new(@values),
                :$nrow,
                :$ncol,
                :$implicit-value
                )
    }

    method transpose() {
        # Make the result arrays
        my CArray[int32] $t-row-ptr = CArray[int32].new(0 xx $!ncol+1);
        my CArray[int32] $t-col-index = CArray[int32].new(0 xx $!col-index.elems);
        my CArray[num64] $t-values = CArray[num64].new(0.0 xx $!values.elems);

        # Call the C routine
        transpose($!row-ptr, $!col-index, $!values, $!nrow, $!ncol, $t-row-ptr, $t-col-index, $t-values);

        # Result
        return CSR.new(
                row-ptr => $t-row-ptr,
                col-index => $t-col-index,
                values => $t-values,
                nrow => $!ncol,
                ncol => $!nrow,
                :$!implicit-value
                );
    }

    method dot-pattern(CSR:D $other) {
        die 'The number of rows of the argument is expected to be equal to the number of columns of the object.'
        unless $!ncol == $other.nrow;

        # Estimate non-zeroes
        #my $dot-nrow = $!row-ptr.Array.rotor(2=>-1).map({ $_.tail > $_.head }).sum;
        my $dot-nrow = $!nrow;
        my $dot-ncol = $other.col-index.Array.unique.elems;
        my $dot-total = $dot-nrow * $dot-ncol;

        # Make the result arrays
        my CArray[int32] $dot-row-ptr = CArray[int32].new(0 xx ($!nrow + 1));
        my CArray[int32] $dot-col-index = CArray[int32].new(0 xx $dot-total);

        # Call the C routine
        dot-pattern(
                $!row-ptr, $!col-index, $other.row-ptr, $other.col-index,
                $!nrow, $!ncol, $other.ncol,
                $dot-row-ptr, $dot-col-index);

        # Result
        return CSR.new(
                row-ptr => $dot-row-ptr,
                col-index => $dot-col-index.head($dot-row-ptr.tail),
                values => (1.0 xx $dot-row-ptr.tail),
                :$!nrow,
                ncol => $other.ncol,
                :$!implicit-value
                );
    }
}


#-----------------------------------------------------------
#`[
our proto sub dot-product(CSR:D $m1, CSR:D $m2 --> CSR:D) is export {*}

multi sub dot-product(CSR:D $m1, CSR:D $m2 --> CSR:D) {
    die 'The number of rows of the argument is expected to be equal to the number of columns of the object.'
    unless $m1.ncol == $m2.nrow;


}
]