name = libSparseMatrixFunctions.dylib

all: /Volumes/Macintosh\ HD/Users/antonov/Raku-Math-SparseMatrix-Native/resources/libraries/libSparseMatrixFunctions.dylib
clean:
	rm /Volumes/Macintosh\ HD/Users/antonov/Raku-Math-SparseMatrix-Native/resources/libraries/libSparseMatrixFunctions.dylib

/Volumes/Macintosh\ HD/Users/antonov/Raku-Math-SparseMatrix-Native/resources/libraries/libSparseMatrixFunctions.dylib: SparseMatrixFunctions.c
ifeq ($(shell uname), Darwin)
	clang -shared  -fno-omit-frame-pointer -fno-optimize-sibling-calls -Werror=vla -O3 -DNDEBUG -Wno-logical-op-parentheses -D_DARWIN_USE_64_BIT_INODE=1  -gdwarf-4 -DMVM_HEAPSNAPSHOT_FORMAT=2 -D_GNU_SOURCE SparseMatrixFunctions.c -o /Volumes/Macintosh\ HD/Users/antonov/Raku-Math-SparseMatrix-Native/resources/libraries/libSparseMatrixFunctions.dylib -framework Accelerate
else
	clang -shared  -fno-omit-frame-pointer -fno-optimize-sibling-calls -Werror=vla -O3 -DNDEBUG -Wno-logical-op-parentheses -D_DARWIN_USE_64_BIT_INODE=1  -gdwarf-4 -DMVM_HEAPSNAPSHOT_FORMAT=2 -D_GNU_SOURCE SparseMatrixFunctions.c -o /Volumes/Macintosh\ HD/Users/antonov/Raku-Math-SparseMatrix-Native/resources/libraries/libSparseMatrixFunctions.dylib
endif