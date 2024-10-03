use v6.d;

unit module Math::SparseMatrix::Native;

use NativeCall;
use NativeHelpers::Array;

my constant $library = %?RESOURCES<libraries/SparseMatrixFunctions>;

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

    sub random_sparse_matrix(CSRStruct is rw, uint32 $nrow, uint32 $ncol, uint32 $nnz, num64 $implicit_value, uint32 $seed --> int32)
            is native($library) {*}

    sub transpose(CSRStruct is rw, CSRStruct --> int32)
            is native($library) {*}

    sub dot_dense_vector(CArray[num64] is rw, CSRStruct, CArray[num64] --> int32 )
            is native($library) {*}

    sub dot_pattern(CSRStruct is rw, CSRStruct, CSRStruct, int32 $nnz --> int32)
            is native($library) {*}

    sub dot_numeric(CSRStruct is rw, CSRStruct, CSRStruct, int32 $nnz --> int32)
            is native($library) {*}

    sub add_scalar_to_sparse_matrix(CSRStruct is rw, CSRStruct $matrix, num64 $scalar, int32 $clone --> int32)
            is native($library) {*}

    sub add_sparse_matrices(CSRStruct is rw, CSRStruct, CSRStruct --> int32)
            is native($library) {*}

    #----------------------------------------------------------------
    method dimensions() { return ($!nrow, $!ncol); }

    #----------------------------------------------------------------
    submethod BUILD(:$values, :$col_index, :$row_ptr, :$nnz is copy = 0,
                    UInt:D :$nrow = 0, UInt:D :$ncol = 0,
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
    multi method random(
            UInt:D $nrow,
            UInt:D $ncol,
            $nnz is copy = Whatever,
            Numeric:D $implicit_value = 0.0,
            :$seed is copy = Whatever --> int32) {
        return self.random(:$nrow, :$ncol, :$nnz, :$implicit_value, :$seed);
    }

    multi method random(UInt:D :$nrow,
                        UInt:D :$ncol,
                        :$nnz is copy = Whatever,
                        Numeric:D :$implicit_value = 0.0,
                        :$seed is copy = Whatever --> int32) {
        if $nnz.isa(Whatever) {
            $nnz = min(1000, $nrow * $ncol * 0.01).Int;
        }
        die 'The agument $nnz is expected to be a non-negative integer or Whatever.'
        unless $nnz ~~ Int:D && $nnz ≥ 0;

        if $seed.isa(Whatever) {
            $seed = now.UInt;
        }
        die 'The agument $seed is expected to be a non-negative integer or Whatever.'
        unless $seed ~~ Int:D && $nnz ≥ 0;


        if $!values || $!col_index || $!row_ptr {
            destroy_sparse_matrix(self);
        }

        my $res = random_sparse_matrix(self, $nrow, $ncol, $nnz, $implicit_value.Num, $seed.UInt);
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
    # Matrix-dense-vector
    multi method dot(@vector) {
        die "If the first argument is a vector, then it is expected to be numeric with length that matches the columns of the sparse matrix."
        unless @vector.elems == $!ncol && (@vector.all ~~ Numeric:D);
        my $vec = CArray[num64].new(@vector);
        my $target = CArray[num64].allocate($!nrow);
        my $res = dot_dense_vector($target, self, $vec);
        return $target.Array;
    }


    #----------------------------------------------------------------
    # Matrix-matrix
    method dot-pattern(CSRStruct $other) {
        my $target = CSRStruct.new(
                values => Nil, col_index => Nil, col_index => Nil,
                :$!nrow, ncol => $other.ncol, nnz => 0, implicit_value => 0e0);
        my $res = dot_pattern($target, self, $other, -1);
        return $target;
    }

    multi method dot(CSRStruct $other) {
        my $target = CSRStruct.new(
                values => Nil, col_index => Nil, col_index => Nil,
                :$!nrow, ncol => $other.ncol, nnz => 0, implicit_value => 0e0);
        my $res = dot_numeric($target, self, $other, -1);
        return $target;
    }

    #----------------------------------------------------------------
    multi method add(Numeric:D $a, Bool:D :$clone = True) {
        if $clone {
            my $target = CSRStruct.new(
                    values => Nil, col_index => Nil, col_index => Nil,
                    :$!nrow, :$!ncol, nnz => 0, :$!implicit_value);
            my $res = add_scalar_to_sparse_matrix($target, self, $a.Num, 1);
            return $target;
        } else {
            my $res = add_scalar_to_sparse_matrix(self, self, $a.Num, 0);
            return self;
        }
    }

    multi method add(CSRStruct $other) {
        my $target = CSRStruct.new();
        my $res = add_sparse_matrices($target, self, $other);
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