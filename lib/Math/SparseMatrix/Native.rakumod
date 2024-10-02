use v6.d;

unit module Math::SparseMatrix::Native;

use NativeCall;
use NativeHelpers::Array;

my constant $library = %?RESOURCES<libraries/SparseMatrixFunctions>;

# void dot_pattern(int *IA, int *JA, int *IB, int *JB, int NP, int NQ, int NR, int *IC, int *JC)
sub dot-pattern(CArray[int32], CArray[int32], CArray[int32], CArray[int32],
                int32, int32, int32,
                CArray[int32], CArray[num64]) is native($library) is symbol('dot_pattern_func') {*}

#sub DotProduct(CArray[num64], CArray[num64], int32 --> num64) is native($library) {*}
#sub DotProductFloat(CArray[num32], CArray[num32], int32 --> num32) is native($library) {*}

sub transpose(CArray[int32], CArray[int32], CArray[num64],
              uint32, uint32,
              CArray[int32], CArray[int32], CArray[num64]) is native($library) {*};

#-----------------------------------------------------------

# Define the CSRStruct struct in Raku
class CSRStruct is repr('CStruct') {
    has CArray[num64] $.values;
    has CArray[int32] $.col_index;
    has CArray[int32] $.row_ptr;
    has uint32 $.nnz;
    has uint32 $.nrow;
    has uint32 $.ncol;
    has num64 $.implicit_value;

    #----------------------------------------------------------------
    # Bind to the create_sparse_matrix function
    sub create_sparse_matrix(CSRStruct is rw, uint32 $nrow, uint32 $ncol, uint32 $nnz, num64 $implicit_value --> int32)
            is native($library) {*}

    # Bind to the destroy_sparse_matrix function
    sub destroy_sparse_matrix(CSRStruct is rw)
            is native($library) {*}

    sub random_sparse_matrix(CSRStruct is rw, uint32 $nrow, uint32 $ncol, uint32 $nnz, num64 $implicit_value --> int32)
            is native($library) {*}

    sub transpose(CSRStruct is rw, CSRStruct --> int32)
            is native($library) {*}

    sub dot_pattern(CSRStruct is rw, CSRStruct, CSRStruct, int32 $nnz --> int32)
            is native($library) {*}

    sub dot_numeric(CSRStruct is rw, CSRStruct, CSRStruct, int32 $nnz --> int32)
            is native($library) {*}

    #----------------------------------------------------------------
    submethod BUILD(:$values, :$col_index, :$row_ptr, :$nnz is copy = 0,
                    UInt:D :$nrow, UInt:D :$ncol,
                    Numeric:D :$implicit_value = 0e0) {

        if $nnz.isa(Whatever) { $nnz = $values.elems }
        die 'The agument $nnz is expected to be a non-negative integer or Whatever.'
        unless $nnz ~~ Int:D && $nnz ≥ 0;

        create_sparse_matrix(self, $nrow, $ncol, $nnz, $implicit_value.Num);

        if $values {
            for ^$nnz -> $i {
                self.values[$i] = $values[$i].Num;
                self.col_index[$i] = $col_index[$i];
            }
            for ^$nrow -> $i {
                self.row_ptr[$i] = $row_ptr[$i].Int
            }
        }
    }

    submethod DESTROY {
        destroy_sparse_matrix(self);
    }


    #----------------------------------------------------------------
    method random(UInt:D $nrow, UInt:D $ncol, $nnz is copy = Whatever, Numeric:D $implicit_value = 0.0 --> int32) {
        if $nnz.isa(Whatever) {
            $nnz = min(1000, $nrow * $ncol * 0.01).Int;
        }
        die 'The agument $nnz is expected to be a non-negative integer or Whatever.'
        unless $nnz ~~ Int:D && $nnz ≥ 0;

        if $!values || $!col_index || $!row_ptr {
            destroy_sparse_matrix(self);
        }

        my $res = random_sparse_matrix(self, $nrow, $ncol, $nnz, $implicit_value.Num);
        return self;
    }

    #----------------------------------------------------------------
    method transpose() {
        my $target = CSRStruct.new(
                values => Nil, col_index => Nil, col_index => Nil,
                nrow => self.ncol, ncol => self.nrow, nnz => self.nnz, :$!implicit_value);
        my $res = transpose($target, self);
        return $target;
    }

    #----------------------------------------------------------------
    method dot-pattern(CSRStruct $other) {
        my $target = CSRStruct.new(
                values => Nil, col_index => Nil, col_index => Nil,
                :$!nrow, ncol => $other.ncol, nnz => 0, implicit_value => 0e0);
        my $res = dot_pattern($target, self, $other, -1);
        return $target;
    }

    method dot(CSRStruct $other) {
        my $target = CSRStruct.new(
                values => Nil, col_index => Nil, col_index => Nil,
                :$!nrow, ncol => $other.ncol, nnz => 0, implicit_value => 0e0);
        my $res = dot_numeric($target, self, $other, -1);
        return $target;
    }

    #----------------------------------------------------------------
    # Method to display/print the sparse matrix
    method print() {
        for ^$!nrow -> $i {
            say "Row ", $i, ":";
            for $!row_ptr[$i] .. $!row_ptr[$i + 1] - 1 -> $j {
                say "  Column ", $!col_index[$j], " = ", $!values[$j];
            }
        }
        say "Implicit value: ", $!implicit_value;
    }
}


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
        my CArray[int32] $t-row-ptr = CArray[int32].new(0 xx $!ncol + 1);
        my CArray[int32] $t-col-index = CArray[int32].new(0 xx $!col-index.elems);
        my CArray[num64] $t-values = CArray[num64].new(0.0 xx $!values.elems);

        # Call the C routine
        transpose($!row-ptr, $!col-index, $!values, $!nrow, $!ncol, $t-row-ptr, $t-col-index, $t-values);

        # Result
        return self.bless(
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