#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Definition of the CSRStruct struct
typedef struct CSRStruct {
    double *values;
    int *col_index;
    int *row_ptr;
    int nnz;
    int nrow;
    int ncol;
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
int create_sparse_matrix(CSRStruct *matrix, int nrow, int ncol, int nnz, double implicit_value) {
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
// Equivalence of two CStructs
//=====================================================================
int eqv_sorted_columns(CSRStruct *matrix1, CSRStruct *matrix2, double tol) {
    if (matrix1->nrow != matrix2->nrow || matrix1->ncol != matrix2->ncol || matrix1->nnz != matrix2->nnz || matrix1->implicit_value != matrix2->implicit_value) {
        return 0;
    }

    for (int i = 0; i < matrix1->nnz; ++i) {
        if (fabs(matrix1->values[i] - matrix2->values[i]) > tol || matrix1->col_index[i] != matrix2->col_index[i]) {
            return 0;
        }
    }

    for (int i = 0; i <= matrix1->nrow; ++i) {
        if (matrix1->row_ptr[i] != matrix2->row_ptr[i]) {
            return 0;
        }
    }

    return 1;
}

int eqv_general(CSRStruct *matrix1, CSRStruct *matrix2, double tol) {
    if (matrix1->nrow != matrix2->nrow || matrix1->ncol != matrix2->ncol || matrix1->nnz != matrix2->nnz) {
        return 0;
    }

    if (fabs(matrix1->implicit_value - matrix2->implicit_value) > tol) {
        return 0;
    }

    for (int i = 0; i < matrix1->nrow; ++i) {
        int start1 = matrix1->row_ptr[i];
        int end1 = matrix1->row_ptr[i + 1];
        int start2 = matrix2->row_ptr[i];
        int end2 = matrix2->row_ptr[i + 1];

        if ((end1 - start1) != (end2 - start2)) {
            return 0;
        }

        for (int j = start1; j < end1; ++j) {
            int found = 0;
            for (int k = start2; k < end2; ++k) {
                if (matrix1->col_index[j] == matrix2->col_index[k] &&
                    fabs(matrix1->values[j] - matrix2->values[k]) <= tol) {
                    found = 1;
                    break;
                }
            }
            if (!found) {
                return 0;
            }
        }
    }

    return 1;
}

//=====================================================================
// Creation from triplets
//=====================================================================
int compare_triplets(const void *a, const void *b) {
    int row_a = ((int *)a)[0];
    int row_b = ((int *)b)[0];
    return row_a - row_b;
}

int create_sparse_matrix_from_triplets(CSRStruct *matrix,
                                        int nrow, int ncol, int nnz,
                                        double implicit_value,
                                        int *rows, int *cols, double *values) {
    if (!matrix || !rows || !cols || !values) return -1;

    int (*triplets)[3] = malloc(nnz * sizeof(*triplets));
    if (!triplets) return -1;

    for (int i = 0; i < nnz; i++) {
        triplets[i][0] = rows[i];
        triplets[i][1] = cols[i];
        triplets[i][2] = i; // Store index for values
    }

    qsort(triplets, nnz, sizeof(*triplets), compare_triplets);

    if (create_sparse_matrix(matrix, nrow, ncol, nnz, implicit_value) != 0) {
        free(triplets);
        return -1;
    }

    matrix->values = (double *)malloc(nnz * sizeof(double));
    matrix->col_index = (int *)malloc(nnz * sizeof(int));
    matrix->row_ptr = (int *)calloc(nrow + 1, sizeof(int));

    if (!matrix->values || !matrix->col_index || !matrix->row_ptr) {
        free(matrix->values);
        free(matrix->col_index);
        free(matrix->row_ptr);
        free(triplets);
        return -1;
    }

    for (int i = 0; i < nnz; i++) {
        int row = triplets[i][0];
        int col = triplets[i][1];
        int val_index = triplets[i][2];
        matrix->values[i] = values[val_index];
        matrix->col_index[i] = col;
        matrix->row_ptr[row + 1]++;
    }

    for (int i = 1; i <= nrow; i++) {
        matrix->row_ptr[i] += matrix->row_ptr[i - 1];
    }

    free(triplets);
    return 0;
}

//=====================================================================
// Random sparse matrix for CStruct
//=====================================================================
int random_sparse_matrix(CSRStruct *matrix, int nrow, int ncol, int nnz, double implicit_value, int seed) {
    srand(seed);

    int *rows = (int *)malloc(nnz * sizeof(int));
    int *cols = (int *)malloc(nnz * sizeof(int));
    double *values = (double *)malloc(nnz * sizeof(double));

    int *used = (int *)calloc(nrow * ncol, sizeof(int));
    int count = 0;

    while (count < nnz) {
        int row = rand() % nrow;
        int col = rand() % ncol;
        if (!used[row * ncol + col]) {
            used[row * ncol + col] = 1;
            rows[count] = row;
            cols[count] = col;
            values[count] = (double)rand() / RAND_MAX; // Random value between 0 and 1
            count++;
        }
    }

    free(used);

    int result = create_sparse_matrix_from_triplets(matrix, nrow, ncol, nnz, implicit_value, rows, cols, values);

    free(rows);
    free(cols);
    free(values);

    return result;
}

//=====================================================================
// Transpose for CStruct
//=====================================================================
int transpose(CSRStruct *target, CSRStruct *matrix) {
    int *IAT = (int *)calloc(matrix->ncol + 1, sizeof(int));
    int *JAT = (int *)malloc(matrix->nnz * sizeof(int));
    double *ANT = (double *)malloc(matrix->nnz * sizeof(double));

    int MH = matrix->ncol + 1;
    int NH = matrix->nrow + 1;

    int IAB = matrix->row_ptr[NH - 1];

    for (int i = 0; i < IAB; ++i) {
        int J = matrix->col_index[i] + 2;
        if (J < MH) {
            IAT[J] += 1;
        }
    }

    IAT[0] = 0;
    IAT[1] = 0;

    if (matrix->ncol != 1) {
        for (int i = 2; i < MH; ++i) {
            IAT[i] += IAT[i - 1];
        }
    }

    for (int i = 0; i < matrix->nrow; ++i) {
        int IAA = matrix->row_ptr[i];
        int IAB = matrix->row_ptr[i + 1];
        if (IAB < IAA) continue;
        for (int jp = IAA; jp < IAB; ++jp) {
            int J = matrix->col_index[jp] + 1;
            int K = IAT[J];
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
// Dot product for CStruct (Matrix-Matrix)
//=====================================================================
int dot_dense_vector(double *target, CSRStruct *matrix, double *vector) {
    for (int i = 0; i < matrix->nrow; i++) {
        target[i] = 0.0;
        int row_start = matrix->row_ptr[i];
        int row_end = matrix->row_ptr[i + 1];
        for (int j = row_start; j < row_end; j++) {
            target[i] += matrix->values[j] * vector[matrix->col_index[j]];
        }
    }
    return 0;
}

//=====================================================================
// Dot product for CStruct (Matrix-Matrix)
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
        int dot_nrow = 0;
        for (int i = 0; i < A->nrow; ++i) {
            if (A->row_ptr[i + 1] > A->row_ptr[i]) {
                dot_nrow++;
            }
        }

        int *unique_cols = (int *)calloc(B->ncol, sizeof(int));
        int dot_ncol = 0;
        for (int i = 0; i < B->nnz; ++i) {
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

    for (int i = 0; i < A->nrow; ++i) {
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
    if (err) { return err; }

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

    for (int i = 0; i < A->nrow; ++i) {
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
    free(pattern.col_index);
    free(pattern.values);
    return 0;
}

//=====================================================================
// Addition-pattern (element-wise)
//=====================================================================
int add_pattern(CSRStruct *result, CSRStruct *matrix, CSRStruct *other) {
    if(matrix->nrow != other->nrow || matrix->ncol != other->ncol) {
        exit(EXIT_FAILURE); // The dimensions of the argument must match the dimensions of the object.
    }

    int *IC = (int*) calloc(matrix->nrow + 1, sizeof(int));
    int *JC = (int*) malloc(matrix->nnz * sizeof(int));
    int *IX = (int*) calloc(matrix->ncol, sizeof(int));
    int IP = 0;

    for(int i = 0; i < matrix->nrow; i++) {
        IC[i] = IP;
        int IAA = matrix->row_ptr[i];
        int IAB = matrix->row_ptr[i + 1] - 1;
        if(IAB >= IAA) {
            for(int jp = IAA; jp <= IAB; jp++) {
                int j = matrix->col_index[jp];
                JC[IP++] = j;
                IX[j] = i + 1;
            }
        }

        int IBA = other->row_ptr[i];
        int IBB = other->row_ptr[i + 1] - 1;
        if(IBB >= IBA) {
            for(int jp = IBA; jp <= IBB; jp++) {
                int j = other->col_index[jp];
                if(IX[j] != i + 1) {
                    JC[IP++] = j;
                }
            }
        }
    }

    IC[matrix->nrow] = IP;

    create_sparse_matrix(result, matrix->nrow, matrix->ncol, IP, 0.0);

    for(int i = 0; i < IP; i++) {
        result->values[i] = 1.0;
    }

    for(int i = 0; i <= matrix->nrow; i++) {
        result->row_ptr[i] = IC[i];
    }

    for(int i = 0; i < IP; i++) {
        result->col_index[i] = JC[i];
    }

    free(IC);
    free(JC);
    free(IX);

    return 0;
}

//=====================================================================
// Addition numeric
//=====================================================================
int add_numeric(CSRStruct *result, CSRStruct *matrix, CSRStruct *other, int op) {
    CSRStruct pattern;
    int err = add_pattern(&pattern, matrix, other);
    if (err) { return err; }

    double *CN = (double*) calloc(pattern.nnz, sizeof(double));
    double *X = (double*) calloc(pattern.ncol, sizeof(double));

    for(int i = 0; i < matrix->nrow; i++) {
        int IH = i + 1;
        int ICA = pattern.row_ptr[i];
        int ICB = pattern.row_ptr[IH] - 1;

        if(ICB < ICA) continue;

        for(int ip = ICA; ip <= ICB; ip++) {
            X[pattern.col_index[ip]] = 0;
        }

        int IAA = matrix->row_ptr[i];
        int IAB = matrix->row_ptr[IH] - 1;

        if(IAB >= IAA) {
            for(int ip = IAA; ip <= IAB; ip++) {
                X[matrix->col_index[ip]] = matrix->values[ip];
            }
        }

        int IBA = other->row_ptr[i];
        int IBB = other->row_ptr[IH] - 1;

        if(IBB >= IBA) {
            for(int ip = IBA; ip <= IBB; ip++) {
                int J = other->col_index[ip];
                X[J] += other->values[ip];
            }
        }

        for(int ip = ICA; ip <= ICB; ip++) {
            CN[ip] = X[pattern.col_index[ip]];
        }
    }

    create_sparse_matrix(result, matrix->nrow, matrix->ncol, pattern.nnz, matrix->implicit_value + other->implicit_value);

    for(int i = 0; i < pattern.nnz; i++) {
        result->values[i] = CN[i];
    }

    for(int i = 0; i <= matrix->nrow; i++) {
        result->row_ptr[i] = pattern.row_ptr[i];
    }

    for(int i = 0; i < pattern.nnz; i++) {
        result->col_index[i] = pattern.col_index[i];
    }

    destroy_sparse_matrix(&pattern);
    free(CN);
    free(X);

    return 0;
}

// TBD
//int multiply_numeric(CSRStruct *result, CSRStruct *matrix, CSRStruct *other) {
//    return op_numeric(result, matrix, other, MULT_OP);
//}

//=====================================================================
// Element-wise generic
//=====================================================================
#define MULT_OP 101
#define ADD_OP 102

int op_scalar_to_sparse_matrix(CSRStruct *result, CSRStruct *matrix, double scalar, int clone, int op) {

    if (clone) {
        result->values = (double*)malloc(matrix->nnz * sizeof(double));
        result->col_index = (int*)malloc(matrix->nnz * sizeof(int));
        result->row_ptr = (int*)malloc((matrix->nrow + 1) * sizeof(int));
        result->nnz = matrix->nnz;
        result->nrow = matrix->nrow;
        result->ncol = matrix->ncol;
        if (op == ADD_OP) {
            result->implicit_value = matrix->implicit_value + scalar;
        } else {
            result->implicit_value = matrix->implicit_value * scalar;
        }

        for (int i = 0; i < matrix->nnz; i++) {
            result->values[i] = matrix->values[i];
            result->col_index[i] = matrix->col_index[i];
            if (op == ADD_OP) {
                result->values[i] += scalar;
            } else {
                result->values[i] *= scalar;
            }
        }
        for (int i = 0; i <= matrix->nrow; i++) {
            result->row_ptr[i] = matrix->row_ptr[i];
        }
    } else {
        if (op == ADD_OP) {
            matrix->implicit_value = matrix->implicit_value + scalar;
            for (int i = 0; i < matrix->nnz; i++) {
                matrix->values[i] += scalar;
            }
        } else {
            matrix->implicit_value = matrix->implicit_value * scalar;
            for (int i = 0; i < matrix->nnz; i++) {
                matrix->values[i] *= scalar;
            }
        }
    }

    return 0;
}


// Note that this routine assumes that the column indexes are sorted per row.
// Hence, in the Raku invoker methods we sort those column indices by calling transpose twice.
int op_sparse_matrices(CSRStruct *result, const CSRStruct *A, const CSRStruct *B, int op) {
    if (A->nrow != B->nrow || A->ncol != B->ncol) return -1;

    int *row_ptr = (int *)calloc(A->nrow + 1, sizeof(int));
    int nnz_estimate = A->nnz + B->nnz;
    double *values = (double *)malloc(nnz_estimate * sizeof(double));
    int *col_index = (int *)malloc(nnz_estimate * sizeof(int));

    int pos = 0;
    for (int i = 0; i < A->nrow; ++i) {
        int a_start = A->row_ptr[i];
        int a_end = A->row_ptr[i + 1];
        int b_start = B->row_ptr[i];
        int b_end = B->row_ptr[i + 1];

        while (a_start < a_end && b_start < b_end) {
            if (A->col_index[a_start] < B->col_index[b_start]) {
                if (op == ADD_OP) {
                    values[pos] = A->values[a_start] + B->implicit_value;
                } else {
                    values[pos] = A->values[a_start] * B->implicit_value;
                }
                col_index[pos] = A->col_index[a_start];
                a_start++;
            } else if (A->col_index[a_start] > B->col_index[b_start]) {
                if (op == ADD_OP) {
                    values[pos] = B->values[b_start] + A->implicit_value;
                } else {
                    values[pos] = B->values[b_start] * A->implicit_value;
                }
                col_index[pos] = B->col_index[b_start];
                b_start++;
            } else {
                if (op == ADD_OP) {
                    values[pos] = A->values[a_start] + B->values[b_start];
                } else {
                    values[pos] = A->values[a_start] * B->values[b_start];
                }
                col_index[pos] = A->col_index[a_start];
                a_start++;
                b_start++;
            }
            pos++;
        }

        while (a_start < a_end) {
            if (op == ADD_OP) {
                values[pos] = A->values[a_start] + B->implicit_value;
            } else {
                values[pos] = A->values[a_start] * B->implicit_value;
            }
            col_index[pos] = A->col_index[a_start];
            a_start++;
            pos++;
        }

        while (b_start < b_end) {
            if (op == ADD_OP) {
                values[pos] = B->values[b_start] + A->implicit_value;
            } else {
                values[pos] = B->values[b_start] * A->implicit_value;
            }
            col_index[pos] = B->col_index[b_start];
            b_start++;
            pos++;
        }

        row_ptr[i + 1] = pos;
    }

    result->values = (double *)realloc(values, pos * sizeof(double));
    result->col_index = (int *)realloc(col_index, pos * sizeof(int));
    result->row_ptr = (int *)realloc(row_ptr, (A->nrow + 1) * sizeof(int));
    result->nnz = pos;
    result->nrow = A->nrow;
    result->ncol = A->ncol;
    if (op == ADD_OP) {
        result->implicit_value = A->implicit_value + B->implicit_value;
    } else {
        result->implicit_value = A->implicit_value * B->implicit_value;
    }

    return 0;
}

//=====================================================================
// Element-wise addition
//=====================================================================
int add_scalar_to_sparse_matrix(CSRStruct *result, CSRStruct *matrix, double scalar, int clone) {
    return op_scalar_to_sparse_matrix(result, matrix, scalar, clone, ADD_OP);
}

int add_sparse_matrices(CSRStruct *result, const CSRStruct *A, const CSRStruct *B) {
    return op_sparse_matrices(result, A, B, ADD_OP);
}

//=====================================================================
// Element-wise multiplication
//=====================================================================
int multiply_scalar_to_sparse_matrix(CSRStruct *result, CSRStruct *matrix, double scalar, int clone) {
    return op_scalar_to_sparse_matrix(result, matrix, scalar, clone, MULT_OP);
}

int multiply_sparse_matrices(CSRStruct *result, const CSRStruct *A, const CSRStruct *B) {
    return op_sparse_matrices(result, A, B, MULT_OP);
}