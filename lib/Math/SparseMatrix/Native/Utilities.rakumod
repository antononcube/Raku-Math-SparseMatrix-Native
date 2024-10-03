unit module Math::SparseMatrix::Native::Utilities;

## Essentially the same as Math::SparseMatrix::Utilities
## but I prefer Math::SparseMatrix::Native to be independent.

#=====================================================================
# Element-wise operation
#=====================================================================
sub elementwise-operation(@matrix1, @matrix2, &op) is export {
    die "Incompatible matrices" unless @matrix1.elems == @matrix2.elems && @matrix1[0].elems == @matrix2[0].elems;
    my @result;
    for @matrix1.kv -> $i, @row {
        for @row.kv -> $j, $elem {
            @result[$i][$j] = op($elem, @matrix2[$i][$j]);
        }
    }
    return @result;
}

#=====================================================================
# Dot Product for dense matrices (arrays)
#=====================================================================
sub dense-dot-product(@matrix1, @matrix2) is export {
    die "Incompatible matrices" unless @matrix1[0].elems == @matrix2.elems;
    my @result;
    for @matrix1.kv -> $i, @row {
        for @matrix2[0].keys -> $j {
            @result[$i][$j] = [+] @row Z* @matrix2.map({ $_[$j] });
        }
    }
    return @result;
}