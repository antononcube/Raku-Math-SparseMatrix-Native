#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


//=====================================================================
// Transpose
//=====================================================================
//
// Input:
// IA, JA, AN          given matrix in RR(C)U.
// N                   number of rows of the matrix.
// M                   number of columns of the matrix.
//
// Output:
// IAT, JAT, ANT   transposed matrix in RR(C)O.
//
void transpose(int *IA, int *JA, double *AN, unsigned int N, unsigned int M, int *IAT, int *JAT, double *ANT) {
    int MH = M + 1;
    int NH = N + 1;
    int IAB, IAA, JAA, JP, J, K;

//    *IAT = (int *)calloc(MH + 1, sizeof(int));
//    *JAT = (int *)malloc((IA[NH - 1]) * sizeof(int));
//    *ANT = (double *)malloc((IA[NH - 1]) * sizeof(double));

    for (int I = 1; I <= M; I++) {
        IAT[I] = 0;
    }

    IAB = IA[NH - 1];
    for (int I = 0; I < IAB; I++) {
        J = JA[I] + 2;
        if (J <= M) {
            IAT[J] = IAT[J] + 1;
        }
    }

    IAT[0] = 0;
    IAT[1] = 0;

    if (M != 1) {
        for (int I = 2; I <= M; I++) {
            IAT[I] = IAT[I] + IAT[I - 1];
        }
    }

    for (int I = 0; I < N; I++) {
        IAA = IA[I];
        IAB = IA[I + 1];
        if (IAB < IAA) {
            continue;
        }
        for (JP = IAA; JP < IAB; JP++) {
            J = JA[JP] + 1;
            K = IAT[J];
            JAT[K] = I;
            ANT[K] = AN[JP];
            IAT[J] = K + 1;
        }
    }
}

//=====================================================================
// Dot product
//=====================================================================
//
//  Input:
//  IA, JA      structure of the first matrix in RR(C)U.
//  IB, JB      structure of the second matrix in RR(C)U.
//
//  NP          number of rows of the first matrix.
//  NQ          number of columns of the first matrix and of rows of the second matrix.
//  NR          number of columns of the second matrix.
//
//  Output:
//  IC, JC      structure of the resulting matrix in RR(C)U.
//
//  Working space:
//  IX          of dimension NR, multiple switch.
//
void dot_pattern(int *IA, int *JA, int *IB, int *JB, int NP, int NQ, int NR, int *IC, int *JC) {
    int *IX = (int *)calloc(NR, sizeof(int));
    int IP = 0;

    for (int I = 0; I < NR; I++) {
        IX[I] = 0;
    }

    for (int I = 0; I < NP; I++) {
        IC[I] = IP;
        int IAA = IA[I];
        int IAB = IA[I + 1] - 1;
        if (IAB < IAA) continue;

        for (int JP = IAA; JP <= IAB; JP++) {
            int J = JA[JP];
            int IBA = IB[J];
            int IBB = IB[J + 1] - 1;
            if (IBB < IBA) continue;

            for (int KP = IBA; KP <= IBB; KP++) {
                int K = JB[KP];
                if (IX[K] == I + 1) continue;
                JC[IP] = K;
                IP++;
                IX[K] = I + 1;
            }
        }
    }

    IC[NP] = IP;
    free(IX);
}

//---------------------------------------------------------------------
//
// Input:
// IA, JA, AN  first given matrix in RR(C)U.
// IB, JB, BN  second given matrix in RR(C)U.
// IC, JC      structure of the resulting matrix in RR(C)U.
// NP          number of rows of the first given matrix and of the resulting matrix.
//
// Output: CN         numerical values of the nonzeros of the resulting matrix.
// Working space: X   expanded accumulator, of dimension equal to the number of columns of the resulting matrix.
//
void dot_numeric(int *IA, int *JA, double *AN, int *IB, int *JB, double *BN, int *IC, int *JC, int NP, double *CN) {
    int NR = IC[NP]; // Number of columns of the resulting matrix
    double *X = (double *)calloc(NR, sizeof(double));

    for (int I = 0; I < NP; I++) {
        int ICA = IC[I];
        int ICB = IC[I + 1] - 1;
        if (ICB < ICA) continue;

        for (int J = ICA; J <= ICB; J++) {
            X[JC[J]] = 0.0;
        }

        int IAA = IA[I];
        int IAB = IA[I + 1] - 1;

        for (int JP = IAA - 1; JP < IAB; JP++) {
            int J = JA[JP];
            double A = AN[JP];
            int IBA = IB[J - 1];
            int IBB = IB[J] - 1;
            if (IBB < IBA) continue;

            for (int KP = IBA - 1; KP < IBB; KP++) {
                int K = JB[KP];
                X[K] += A * BN[KP];
            }
        }

        for (int J = ICA; J <= ICB; J++) {
            CN[J] = X[JC[J]];
        }
    }

    free(X);
}

