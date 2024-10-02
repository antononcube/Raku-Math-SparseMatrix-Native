#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// Definition of the CSRStruct struct
typedef struct CSRStruct {
    double *values;
    int *col_index;
    int *row_ptr;
    unsigned int nnz;
    unsigned int nrow;
    unsigned int ncol;
    double implicit_value;
} CSRStruct;

/**
 * @brief Creates and initializes a CSRStruct in CSR (Compressed Sparse Row) format.
 *
 * @param nrow            Number of rows in the sparse matrix.
 * @param ncol            Number of columns in the sparse matrix.
 * @param nnz             Number of non-zero elements in the sparse matrix.
 * @param implicit_value  The implicit value for elements not explicitly stored.
 * @param matrix          Pointer to the CSRStruct struct to be initialized.
 *
 * @return void
 *
 * @note If memory allocation fails for any of the arrays, the corresponding pointers are set to NULL,
 *       and the nnz, nrow, and ncol are set to 0.
 */
int create_sparse_matrix(CSRStruct *matrix, unsigned int nrow, unsigned int ncol, unsigned int nnz, double implicit_value) {
    if (matrix == NULL) {
        return 1;
    }
    
    // Initialize the struct fields
    matrix->nrow = nrow;
    matrix->ncol = ncol;
    matrix->nnz = nnz;
    matrix->implicit_value = implicit_value;

    // Allocate memory for the values array
    if (nnz > 0) {
        matrix->values = (double *)malloc(nnz * sizeof(double));
        if (matrix->values == NULL) {
            // Memory allocation failed
            matrix->nnz = 0;
            matrix->nrow = 0;
            matrix->ncol = 0;
            return 2;
        }
        // Initialize values to zero (or any other desired default)
        memset(matrix->values, 0, nnz * sizeof(double));
    } else {
        matrix->values = NULL;
    }

    // Allocate memory for the column indices array
    if (nnz > 0) {
        matrix->col_index = (int *)malloc(nnz * sizeof(int));
        if (matrix->col_index == NULL) {
            // Memory allocation failed
            free(matrix->values);
            matrix->values = NULL;
            matrix->nnz = 0;
            matrix->nrow = 0;
            matrix->ncol = 0;
            return 2;
        }
        // Initialize column indices to zero
        memset(matrix->col_index, 0, nnz * sizeof(int));
    } else {
        matrix->col_index = NULL;
    }

    // Allocate memory for the row pointers array
    if (nrow > 0) {
        matrix->row_ptr = (int *)malloc((nrow + 1) * sizeof(int));
        if (matrix->row_ptr == NULL) {
            // Memory allocation failed
            free(matrix->values);
            free(matrix->col_index);
            matrix->values = NULL;
            matrix->col_index = NULL;
            matrix->nnz = 0;
            matrix->nrow = 0;
            matrix->ncol = 0;
            return 2;
        }
        // Initialize row pointers to zero
        memset(matrix->row_ptr, 0, (nrow + 1) * sizeof(int));
    } else {
        matrix->row_ptr = NULL;
    }

    return 0;
}

/**
 * @brief Frees the memory allocated for a CSRStruct struct.
 *
 * @param matrix Pointer to the CSRStruct struct to be freed.
 */
void destroy_sparse_matrix(CSRStruct *matrix) {
    if (matrix == NULL) {
        return;
    }

    // Free the values array
    if (matrix->values != NULL) {
        free(matrix->values);
        matrix->values = NULL;
    }

    // Free the column indices array
    if (matrix->col_index != NULL) {
        free(matrix->col_index);
        matrix->col_index = NULL;
    }

    // Free the row pointers array
    if (matrix->row_ptr != NULL) {
        free(matrix->row_ptr);
        matrix->row_ptr = NULL;
    }

    // Reset other fields
    matrix->nnz = 0;
    matrix->nrow = 0;
    matrix->ncol = 0;
    matrix->implicit_value = 0.0;
}


//=====================================================================
// Random sparse matrix for CStruct
//=====================================================================
int random_sparse_matrix(CSRStruct *matrix, unsigned int nrow, unsigned int ncol, unsigned int nnz, double implicit_value, unsigned int seed) {
    if (matrix == NULL || nrow == 0 || ncol == 0 || nnz == 0 || nnz > nrow * ncol) {
        return 1;
    }

    matrix->values = (double *)malloc(nnz * sizeof(double));
    matrix->col_index = (int *)malloc(nnz * sizeof(int));
    matrix->row_ptr = (int *)calloc(nrow + 1, sizeof(int));
    matrix->nnz = nnz;
    matrix->nrow = nrow;
    matrix->ncol = ncol;
    matrix->implicit_value = implicit_value;

    if (matrix->values == NULL || matrix->col_index == NULL || matrix->row_ptr == NULL) {
        free(matrix->values);
        free(matrix->col_index);
        free(matrix->row_ptr);
        return 1;
    }

    //srand(time(NULL));
    srand(seed);

    for (unsigned int i = 0; i < nnz; ++i) {
        matrix->values[i] = (double)rand() / RAND_MAX;
        matrix->col_index[i] = rand() % ncol;
    }

    for (unsigned int i = 0; i < nnz; ++i) {
        unsigned int row = rand() % nrow;
        ++matrix->row_ptr[row + 1];
    }

    for (unsigned int i = 1; i <= nrow; ++i) {
        matrix->row_ptr[i] += matrix->row_ptr[i - 1];
    }

    return 0;
}

//=====================================================================
// Transpose for CStruct
//=====================================================================
int transpose(CSRStruct *target, CSRStruct *matrix) {
    int *IAT = (int *)calloc(matrix->ncol + 1, sizeof(int));
    int *JAT = (int *)malloc(matrix->nnz * sizeof(int));
    double *ANT = (double *)malloc(matrix->nnz * sizeof(double));

    unsigned int MH = matrix->ncol + 1;
    unsigned int NH = matrix->nrow + 1;

    unsigned int IAB = matrix->row_ptr[NH - 1];

    for (unsigned int i = 0; i < IAB; ++i) {
        unsigned int J = matrix->col_index[i] + 2;
        if (J < MH) {
            IAT[J] += 1;
        }
    }

    IAT[0] = 0;
    IAT[1] = 0;

    if (matrix->ncol != 1) {
        for (unsigned int i = 2; i < MH; ++i) {
            IAT[i] += IAT[i - 1];
        }
    }

    for (unsigned int i = 0; i < matrix->nrow; ++i) {
        unsigned int IAA = matrix->row_ptr[i];
        unsigned int IAB = matrix->row_ptr[i + 1];
        if (IAB < IAA) continue;
        for (unsigned int jp = IAA; jp < IAB; ++jp) {
            unsigned int J = matrix->col_index[jp] + 1;
            unsigned int K = IAT[J];
            JAT[K] = i;
            ANT[K] = matrix->values[jp];
            IAT[J] = K + 1;
        }
    }

    target->values = ANT;
    target->col_index = JAT;
    target->row_ptr = IAT;
    target->nnz = matrix->nnz;
    target->nrow = matrix->ncol;
    target->ncol = matrix->nrow;
    target->implicit_value = matrix->implicit_value;

    return 0;
}

//=====================================================================
// Dot product for CStruct
//=====================================================================

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

int dot_pattern(CSRStruct *result, const CSRStruct *A, const CSRStruct *B, int nnz) {
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
    // This should be refactored to use create_sparse_matrix
    //create_sparse_matrix(result, A->nrow, B->ncol, 0, A->implicit_value);
    result->values = (double *)malloc(IP * sizeof(double));
    for (int i = 0; i < IP; ++i) {
        result->values[i] = 1.0;
    }

    // result->col_index = JC;
    result->col_index = (int *)calloc(IP, sizeof(int));
    for (int i = 0; i < IP; ++i) {
        result->col_index[i] = JC[i];
    }

    result->row_ptr = IC;
    result->nnz = IP;
    result->nrow = A->nrow;
    result->ncol = B->ncol;
    result->implicit_value = 0.0;

    free(IX);
    free(JC);
    return 0;
}

int dot_numeric(CSRStruct *result, const CSRStruct *A, const CSRStruct *B, int nnz) {
    if (A->ncol != B->nrow) {
        fprintf(stderr, "The number of rows of the argument is expected to be equal to the number of columns of the object.\n");
        exit(EXIT_FAILURE);
    }

    CSRStruct pattern;
    int err = dot_pattern(&pattern, A, B, nnz);
    int *IC = pattern.row_ptr;
    int *JC = (int *)calloc(pattern.nnz, sizeof(int));
    for (int i = 0; i < pattern.nnz; ++i) {
        JC[i] = pattern.col_index[i];
    }
    int *IB = B->row_ptr;
    int *JB = B->col_index;
    double *BN = B->values;
    double *X = (double *)calloc(B->ncol, sizeof(double));
    double *result_values = (double *)calloc(pattern.nnz, sizeof(double));
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

    result->values = result_values;
    result->col_index = JC;
    result->row_ptr = IC;
    result->nnz = pattern.row_ptr[pattern.nrow];
    result->nrow = A->nrow;
    result->ncol = B->ncol;
    result->implicit_value = 0.0;

    free(X);
    free(pattern.row_ptr);
    free(pattern.col_index);
    free(pattern.values);
    return 0;
}