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
    has num64 $.implicit_value is rw;

    #----------------------------------------------------------------
    # Bind to the create_sparse_matrix function
    sub create_sparse_matrix(CSRStruct is rw, uint32 $nrow, uint32 $ncol, uint32 $nnz, num64 $implicit_value --> int32)
            is native($library) {*}

    # Bind to the destroy_sparse_matrix function
    sub destroy_sparse_matrix(CSRStruct is rw)
            is native($library) {*}

    sub random_sparse_matrix(CSRStruct is rw, uint32 $nrow, uint32 $ncol, uint32 $nnz, num64 $implicit_value,
                             uint32 $seed --> int32)
            is native($library) {*}

    sub eqv_sorted_columns(CSRStruct, CSRStruct, num64 --> int32)
            is native($library) {*}

    sub eqv_general(CSRStruct, CSRStruct, num64 --> int32)
            is native($library) {*}

    sub transpose(CSRStruct is rw, CSRStruct --> int32)
            is native($library) {*}

    sub dot_dense_vector(CArray[num64] is rw, CSRStruct, CArray[num64] --> int32)
            is native($library) {*}

    sub dot_pattern(CSRStruct is rw, CSRStruct, CSRStruct, int32 $nnz --> int32)
            is native($library) {*}

    sub dot_numeric(CSRStruct is rw, CSRStruct, CSRStruct, int32 $nnz --> int32)
            is native($library) {*}

    sub add_scalar_to_sparse_matrix(CSRStruct is rw, CSRStruct $matrix, num64 $scalar, int32 $clone --> int32)
            is native($library) {*}

    sub add_sparse_matrices(CSRStruct is rw, CSRStruct, CSRStruct --> int32)
            is native($library) {*}

    sub multiply_scalar_to_sparse_matrix(CSRStruct is rw, CSRStruct $matrix, num64 $scalar, int32 $clone --> int32)
            is native($library) {*}

    sub multiply_sparse_matrices(CSRStruct is rw, CSRStruct, CSRStruct --> int32)
            is native($library) {*}

    #----------------------------------------------------------------
    method dimensions() {
        return ($!nrow, $!ncol);
    }

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
    multi method new(:$values!, :$col_index, :$row_ptr, :$nnz is copy = 0,
                     UInt:D :$nrow = 0, UInt:D :$ncol = 0,
                     Numeric:D :$implicit_value = 0e0) {
        self.bless(:$values, :$col_index, :$row_ptr, :$nrow, :$ncol, :$nnz, :$implicit_value);
    }

    multi method new(:dense_matrix(:@dense-matrix)! where @dense-matrix ~~ List:D && @dense-matrix.all ~~ List:D,
                     :$nrow is copy = @dense-matrix.elems,
                     :$ncol is copy = @dense-matrix>>.elems.max,
                     Numeric:D :implicit_value(:$implicit-value) = 0) {
        my @values;
        my @col-index;
        my @row-ptr = 0;
        my $nnz = 0;

        for @dense-matrix.kv -> $row, @cols {
            for @cols.kv -> $col, $val {
                if $val != $implicit-value && $row < $nrow && $col < $ncol {
                    @values.push: $val;
                    @col-index.push: $col;
                    $nnz++;
                }
            }
            @row-ptr.push: @values.elems;
        }

        # Empty trailing rows
        if $nrow > @dense-matrix.elems {
            @row-ptr.append(@row-ptr.tail xx ($nrow - @dense-matrix.elems))
        }

        # Consistency with CSR format
        @row-ptr.append(@row-ptr.tail);

        # Result
        self.bless(
                :@values, col_index => @col-index, row_ptr => @row-ptr,
                :$nrow, :$ncol, :$nnz,
                implicit_value => $implicit-value);
    }

    #----------------------------------------------------------------
    # Clone
    #----------------------------------------------------------------
    method clone() {
        self.new(values => $!values.clone,
                col_index => $!col_index.clone,
                row_ptr => $!row_ptr.clone,
                :$!nrow, :$!ncol, :$!nnz, :$!implicit_value);
    }

    #----------------------------------------------------------------
    # Dense array
    #----------------------------------------------------------------
    #| (Dense) array of arrays representation.
    #| C<:$implicit-value> -- Implicit value to use.
    method Array(:i(:iv(:$implicit-value)) is copy = Whatever) {
        if $implicit-value.isa(Whatever) { $implicit-value = $!implicit_value }
        my @result;
        for ^$!nrow -> $i {
            my @row = ($implicit-value xx $!ncol);
            for $!row_ptr[$i] ..^ $!row_ptr[$i + 1] -> $j {
                @row[$!col_index[$j]] = $!values[$j];
            }
            @result.push(@row);
        }
        return @result;
    }

    #----------------------------------------------------------------
    method eqv(CSRStruct $other, :$method = Whatever, Numeric:D :$tol = 1e-14 --> Bool:D) {
        return do given $method {

            when $_.isa(Whatever) || ($_ ~~ Str:D) && $_.lc ∈ <generic general> {
                eqv_general(self, $other, $tol).so
            }

            when ($_ ~~ Str:D) && $_.lc ∈ <tr transpose transposing> {
                my $m1 = self.transpose.transpose;
                my $m2 = $other.transpose.transpose;
                eqv_sorted_columns($m1, $m2, $tol).so
            }

            default {
                die 'The value of the argument $method is expected to be "generic", "transposing", or Whatever.'
            }
        }
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
        die "If the first argument is a vector, then it is expected to be numeric, with length that matches the columns of the sparse matrix."
        unless @vector.elems == $!ncol && (@vector.all ~~ Numeric:D);
        my $vec = CArray[num64].new(@vector);
        my $target = CArray[num64].allocate($!nrow);
        my $res = dot_dense_vector($target, self, $vec);
        #return self.new(dense-matrix => $target.Array.map({ [$_,] }), :$!nrow, ncol => 1, nnz => $!nrow, implicit_value => 0);
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
    multi method multiply(Numeric:D $a, Bool:D :$clone = True) {
        if $clone {
            my $target = CSRStruct.new(
                    values => Nil, col_index => Nil, col_index => Nil,
                    :$!nrow, :$!ncol, nnz => 0, :$!implicit_value);
            my $res = multiply_scalar_to_sparse_matrix($target, self, $a.Num, 1);
            return $target;
        } else {
            my $res = multiply_scalar_to_sparse_matrix(self, self, $a.Num, 0);
            return self;
        }
    }

    multi method multiply(CSRStruct $other) {
        my $target = CSRStruct.new();
        my $res = multiply_sparse_matrices($target, self, $other);
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