//
//  homography.c
//  LUdecomposition
//
//  Created by Emanuel Brici on 12/2/15.
//  Copyright © 2015 Emanuel Brici. All rights reserved.


#include "LUdecomp.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


// Prints out the matrix in row major order
void printMatrix(double **H) {
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            printf("%lf ", H[i][j]);
        }
        printf("\n");
    }
    
}

void homography(double **A, double *B, int N) {
    double X[8];
    double **H = createMatrix(3);
    double **ATA = createMatrix(8);
    double ATB[8];
    int c = 0;
    
    //when N is equal to 4 this method is used
    if(N == 4){
        LUdecomp* ADecomp = LUdecompose(2 * N, A);
        LUsolve(ADecomp, B, X);
        LUdestroy(ADecomp);
        //create matrix h used for printing homography
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++) {
                if(i == 2 && j == 2) {
                    H[i][j] = 1;
                }
                else {
                    H[i][j] = X[c];
                    c++;
                }
            }
        }
        
    //When N is greater then 4 there is an over determined  system and this method is used instead
    } else {
        double sum = 0.0;
        for(int i = 0; i < 8; i++){
            for(int j = i; j < 8; j++){
                sum = 0.0;
                for(int k = 0; k < 2*N; k++){
                    sum += A[k][i]*A[k][j];
                }
                ATA[i][j] = sum;
            }
        }
        for(int i = 0; i < 8; i++)
            for(int j = 0; j < i; j++)
                ATA[i][j] = ATA[j][i];
        
        sum = 0.0;
        for(int i = 0; i < 8; i++){
            sum = 0.0;
            for(int j = 0; j < 2*N; j++)
                sum += A[j][i] * B[j];
            ATB[i] = sum;
        }
        
        LUdecomp* ATDecomp = LUdecompose(8, ATA);
        LUsolve(ATDecomp, ATB, X);
        LUdestroy(ATDecomp);
        //create matrix h used for printing homography
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++) {
                if(i == 2 && j == 2) {
                    H[i][j] = 1;
                }
                else {
                    H[i][j] = X[c];
                    c++;
                }
            }
        }
        
    }
    printMatrix(H);
}

void matrixSetUp(double **A, double*B, int N){
    double startX[N];
    double startY[N];
    double finishX[N];
    double finishY[N];
    
    //scan in x and y for the first matrix
    for(int i = 0; i < N; i++) {
        scanf("%lf %lf", &startX[i], &startY[i]);
    }
    //scan in x and y for the second matrix
    for(int i = 0; i < N; i++) {
        scanf("%lf %lf", &finishX[i], &finishY[i]);
    }
    
    //fill in matrix how was shown in equation 7 from documentation
    for(int i = 0, j = 0; i < N *2 ; i++, j++) {
        A[i][0] = startX[j];
        A[i][1] = startY[j];
        A[i][2] = 1;
        A[i][3] = 0;
        A[i][4] = 0;
        A[i][5] = 0;
        A[i][6] = (-1) * startX[j] * finishX[j];
        A[i][7] = (-1) * startY[j] * finishX[j];
        B[i] = finishX[j];
        i = i + 1;
        A[i][0] = 0;
        A[i][1] = 0;
        A[i][2] = 0;
        A[i][3] = startX[j];
        A[i][4] = startY[j];
        A[i][5] = 1;
        A[i][6] = (-1) * startX[j] * finishY[j];
        A[i][7] = (-1) * startY[j] * finishY[j];
        B[i] = finishY[j];
    }
    
}

int main(){
    int N;
     //N is number of source and target points
    scanf("%d", &N);
    double **A = createMatrix(N*2);
    double B[2 * N];
    matrixSetUp(A, B, N);
    homography(A,B,N);
    return 0;
}

    
    

