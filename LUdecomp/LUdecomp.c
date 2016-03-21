//
//  LUdecomp.c
//  LUdecomposition
//
//  Created by Emanuel Brici on 11/9/15.
//  Copyright © 2015 Emanuel Brici. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LUdecomp.h"

double **createMatrix(int N) {
    double **M = (double **) malloc(N*sizeof(double*));
    for (int i = 0; i < N; i++)
        M[i] = (double*) malloc(N*sizeof(double));
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            M[i][j] = (i == j) ? 1.0 : 0.0;
    return M;
}

void LUdestroy(LUdecomp* M) {
    for (int i = 0; i < M->N; i++)
        free(M->LU[i]);
    free(M->LU);
}

LUdecomp *LUdecompose(int N, const double **A) {
    LUdecomp *LU = (LUdecomp*) malloc(sizeof(LUdecomp));
    LU->N = N;
    LU->LU = createMatrix(N);
    LU->mutate = (short *) malloc(N*sizeof(short));       //Initialize row permutation array
    LU->d = +1;                                           //Initialize row swap parity value.
    
    // Clone A into LU
    double **A_ = LU->LU;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            A_[i][j] = A[i][j];
        }
    
    for (int i = 0; i < N; i++)
        LU->mutate[i] = (short) i;
    
    //replace A with LU column by column
    for (int j = 0; j <= N - 1 ; j++) {
        for (int i = 0; i <= j; i++) {                  //Compute aij ← βij on and above diagonal.
            double sum = 0.0;
            for (int k = 0; k <= i -1; k++) {
                sum += A_[i][k]*A_[k][j];
            }
            A_[i][j] = A_[i][j] - sum;
        }
        double p = fabs(A_[j][j]);                      //p = initial pivot value
        int n = j;                                      //n = initial pivot row
        
        for (int i = j + 1 ; i <= N - 1; i++) {         //Compute aij ← αij below diagonal.
            double sum = 0.0;
            for (int k = 0; k <= j - 1; k++) {
                sum += A_[i][k] * A_[k][j];
            }
            A_[i][j] = A_[i][j] - sum;
            if (fabs(A_[i][j]) > p) {                   //If better pivot found
                p = fabs(A_[i][j]);                     //record new pivot.
                n = i;
            }
            
        }
        if (p == 0) {                                   //Singular matrix! If p ≈ 0 we may have problems.
            printf("P is equal to 0");
            exit (1);
        }
        if (n != j) {                                   // If best pivot off diagonal
            short temp = LU->mutate[n];
            double *tmp = (double*) malloc(N*sizeof(double));
            for (int i = 0; i < N; i++) {               //swap rows n and j of A
                tmp[i] = A_[n][i];
                A_[n][i] = A_[j][i];
                A_[j][i] = tmp[i];
            }
            temp = LU->mutate[n];                       //swap mutate[n] and mutate[j]
            LU->mutate[n] = LU->mutate[j];
            LU->mutate[j] = temp;
            LU->d = -(LU->d);
        }
        for (int i = j + 1 ; i <= N-1; i++) {           // perform divisions below diagonal
            A_[i][j] = (A_[i][j])/(A_[j][j]);
        }
    }
    return LU;                                          //return decomposed matrix
}

void LUsolve(LUdecomp *decomp, const double *b, double *x){
    double **deCompMatrix = decomp->LU;
    int n = decomp->N;
    double *y = (double*) malloc(n * sizeof(double));
    y[0] = b[decomp->mutate[0]];
    for (int i = 1; i < n; i++) {                       //get zeros below diagonal
        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += deCompMatrix[i][j] * y[j];
        }
        y[i] = (b[decomp->mutate[i]]) - sum;
    }
    x[n - 1] = (y[n - 1]) / (deCompMatrix[n - 1][n - 1]);  //backsolve
    for (int i = n - 2; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += (deCompMatrix[i][j]) * (x[j]);
        }
        x[i] = (1 / deCompMatrix[i][i]) * (y[i] - sum);  //solution to matrix is in x[i]
    }
}










