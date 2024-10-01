#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct {
    double *values;
    int *col_index;
    int *row_ptr;
    unsigned int nnz;
    unsigned int nrow;
    unsigned int ncol;
    double implicit_value;
} SparseMatrixCSR;

SparseMatrixCSR generate_sparse_matrix(unsigned int nrow, unsigned int ncol, unsigned int nnz, double implicit_value) {
    SparseMatrixCSR matrix;
    matrix.nrow = nrow;
    matrix.ncol = ncol;
    matrix.nnz = nnz;
    matrix.implicit_value = implicit_value;
    matrix.values = (double *)malloc(nnz * sizeof(double));
    matrix.col_index = (int *)malloc(nnz * sizeof(int));
    matrix.row_ptr = (int *)malloc((nrow + 1) * sizeof(int));

    srand(time(NULL));
    for (unsigned int i = 0; i < nnz; i++) {
        matrix.values[i] = (double)rand() / RAND_MAX;
        matrix.col_index[i] = rand() % ncol;
    }

    matrix.row_ptr[0] = 0;
    for (unsigned int i = 1; i <= nrow; i++) {
        matrix.row_ptr[i] = matrix.row_ptr[i - 1] + (rand() % (nnz / nrow));
        if (matrix.row_ptr[i] > nnz) {
            matrix.row_ptr[i] = nnz;
        }
    }

    return matrix;
}

void free_sparse_matrix(SparseMatrixCSR *matrix) {
    free(matrix->values);
    free(matrix->col_index);
    free(matrix->row_ptr);
}

int dot_nrow(const int *row_ptr, int n) {
    int dot_nrow = 0;
    for (int i = 0; i < n - 1; i++) {
        if (row_ptr[i + 1] > row_ptr[i]) {
            dot_nrow++;
        }
    }
    return dot_nrow;
}

int dot_ncol(const int *col_index, int nnz) {
    int *unique_col = (int *)malloc(nnz * sizeof(int));
    int unique_count = 0;

    for (int i = 0; i < nnz; i++) {
        int col = col_index[i];
        int is_unique = 1;
        for (int j = 0; j < unique_count; j++) {
            if (unique_col[j] == col) {
                is_unique = 0;
                break;
            }
        }
        if (is_unique) {
            unique_col[unique_count++] = col;
        }
    }

    free(unique_col);
    return unique_count;
}

SparseMatrixCSR dot_pattern(const SparseMatrixCSR *A, const SparseMatrixCSR *B, int nnz) {
    if (A->ncol != B->nrow) {
        fprintf(stderr, "The number of rows of the argument is expected to be equal to the number of columns of the object.\n");
        exit(EXIT_FAILURE);
    }

    if (nnz < 0) {
        unsigned int dot_nrow = 0;
        for (unsigned int i = 0; i < A->nrow; ++i) {
            if (A->row_ptr[i + 1] > A->row_ptr[i]) {
                dot_nrow++;
            }
        }

        unsigned int *unique_cols = (unsigned int *)calloc(B->ncol, sizeof(unsigned int));
        unsigned int dot_ncol = 0;
        for (unsigned int i = 0; i < B->nnz; ++i) {
            if (!unique_cols[B->col_index[i]]) {
                unique_cols[B->col_index[i]] = 1;
                dot_ncol++;
            }
        }
        free(unique_cols);

        nnz = dot_nrow * dot_ncol;
    }

    if (nnz <= 0) {
        fprintf(stderr, "The argument nnz is expected a positive integer or Whatever.\n");
        exit(EXIT_FAILURE);
    }

    int *IC = (int *)calloc(A->nrow + 1, sizeof(int));
    int *JC = (int *)malloc(nnz * sizeof(int));
    int *IX = (int *)calloc(B->ncol, sizeof(int));
    int IP = 0;

    for (unsigned int i = 0; i < A->nrow; ++i) {
        IC[i] = IP;
        int IAA = A->row_ptr[i];
        int IAB = A->row_ptr[i + 1] - 1;
        if (IAB >= IAA) {
            for (int jp = IAA; jp <= IAB; ++jp) {
                int j = A->col_index[jp];
                int IBA = B->row_ptr[j];
                int IBB = B->row_ptr[j + 1] - 1;
                if (IBB >= IBA) {
                    for (int kp = IBA; kp <= IBB; ++kp) {
                        int k = B->col_index[kp];
                        if (IX[k] != i + 1) {
                            JC[IP++] = k;
                            IX[k] = i + 1;
                        }
                    }
                }
            }
        }
    }

    IC[A->nrow] = IP;

    SparseMatrixCSR result;
    result.values = (double *)malloc(IP * sizeof(double));
    for (int i = 0; i < IP; ++i) {
        result.values[i] = 1.0;
    }
    result.col_index = JC;
    result.row_ptr = IC;
    result.nnz = IP;
    result.nrow = A->nrow;
    result.ncol = B->ncol;
    result.implicit_value = 0.0;

    free(IX);
    return result;
}

SparseMatrixCSR dot_numeric(const SparseMatrixCSR *A, const SparseMatrixCSR *B, int nnz) {
    if (A->ncol != B->nrow) {
        fprintf(stderr, "The number of rows of the argument is expected to be equal to the number of columns of the object.\n");
        exit(EXIT_FAILURE);
    }

    SparseMatrixCSR pattern = dot_pattern(A, B, nnz);
    int *IC = pattern.row_ptr;
    int *JC = pattern.col_index;
    int *IB = B->row_ptr;
    int *JB = B->col_index;
    double *BN = B->values;
    double *X = (double *)calloc(B->ncol, sizeof(double));
    double *result_values = (double *)calloc(pattern.row_ptr[pattern.nnz], sizeof(double));
    int IP = 0;

    for (unsigned int i = 0; i < A->nrow; ++i) {
        int ICA = IC[i];
        int ICB = IC[i + 1];

        if (ICB <= ICA) continue;

        for (int j = ICA; j < ICB; ++j) {
            X[JC[j]] = 0;
        }

        int IAA = A->row_ptr[i];
        int IAB = A->row_ptr[i + 1];
        for (int jp = IAA; jp < IAB; ++jp) {
            int j = A->col_index[jp];
            double a = A->values[jp];
            int IBA = IB[j];
            int IBB = IB[j + 1];

            if (IBB <= IBA) continue;

            for (int kp = IBA; kp < IBB; ++kp) {
                int k = JB[kp];
                X[k] += a * BN[kp];
            }
        }

        for (int j = ICA; j < ICB; ++j) {
            result_values[IP++] = X[JC[j]];
        }
    }

    SparseMatrixCSR result;
    result.values = result_values;
    result.col_index = JC;
    result.row_ptr = IC;
    result.nnz = pattern.row_ptr[pattern.nrow];
    result.nrow = A->nrow;
    result.ncol = B->ncol;
    result.implicit_value = 0.0;

    free(X);
    return result;
}

int main() {
    unsigned int nrow = 1000, ncol = 1000, nnz = 10000;
    SparseMatrixCSR A, B, C;

    clock_t start, end;

    start = clock();
    A = generate_sparse_matrix(nrow, ncol, nnz, 0.0);
    B = generate_sparse_matrix(ncol, nrow, nnz, 0.0);
    end = clock();
    printf("Matrix density A: %f\n", (double)(A.nnz) / (double)(A.nrow * A.ncol));
    printf("Matrix density B: %f\n", (double)(B.nnz) / (double)(B.nrow * B.ncol));
    printf("Matrix generation time: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

    start = clock();
    //C = dot_pattern(&A, &B, -1);
    C = dot_numeric(&A, &B, -1);
    end = clock();
    printf("Dot product time: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

    free_sparse_matrix(&A);
    free_sparse_matrix(&B);
    free_sparse_matrix(&C);

    return 0;
}