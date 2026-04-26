
## SVD implementation

The Singular Value Decomposition (SVD) implementation was done with OpenAI's Codex, model "gpt-5.5".
(Several interactions were made in order to get good enough sparse SVD results.)

### First

> Examine the file "./src/SparseMatrixFunctions.c" and implement the routine `svd` the does thin Singular Value Decomposition (SVD) over sparse matrices. 
> In the input is an object CSRStruct, the results matrices are "u", "s", "v". 
> Make sure:
> - You use only data structures and operations implemented in this package; do not use other libraries
> - The sparse SVD implementation is robust and standard
> - The routine `svd`  has an integer argument `k` and returns the decomposition that goes with the `k` largest singular values

### After [profiling](../examples/svd-profiling.raku)

> Your implementation is slow -- hundreds of times slower than expected. It seems you implemented a dense matrix SVD. Maybe, because I said "standard". 
> Implement a "standard sparse matrix SVD".

