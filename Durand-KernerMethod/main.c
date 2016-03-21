//
//  main.c
//  dk
//
//  Created by Emanuel Brici on 10/4/15.
//  Copyright © 2015 Emanuel Brici. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <float.h>

#define MAXN 25 // max degree
#define NUMITER 50 //number of iterations

int main(void) {
    int n = 0;  // polynomial degree
    complex double c[MAXN+1];
    double a, b;
    while (scanf("%lf %lf", &a, &b) == 2)
        c[n++] = a + b*I;
    c[n] = 1;  // leading coeff implied
    
    //find max coefficient
    double max = 0;
    for (int j = 0; j < n; j++) {
        double v = cabs(c[j]);                  //take absolute value of c[j]
        if(v > max){                            //set v to maximum number
            max = v;
        }
    }
    //find radius
    double R = 1 + max;
    double Dtheta = (2*M_PI)/n;                 //divide unit circle by n number of coefficients
    double theta = 0;
    complex double z[n];                        // holds guesses that have been
    for(int j = 0; j<n; j++) {                  //placed around the unit circle
        z[j] = (cos(theta)+ I*sin(theta))*R;    //guess found to be placed in list
        theta += Dtheta;
    }
    
    printf("iter 1 \n");                        //print out first guess
    for (int j = 0; j < n; j++) {
        printf("z[%i] = %0.10f + %0.10f i\n", j, creal(z[j]), cimag(z[j]));
    }
    complex double deltaz[n];                   //hold array of delta z
    for (int k = 1; k < 50; k++) {
        complex double zmax = 0;                //set zmax to 0
        for(int j = 0; j < n; j++) {
            complex double fofx =  1;           //find what f(x) is equal to
            for (int i = n-1; i >= 0; i--) {    //using horner's method to
                fofx = c[i] + fofx * z[j];      //break up the polynomial
            }
            complex double q = 1;               //Q is needed to up date
            for (int i = 0; i < n; i++) {       //the equation that iteratively
                if (i != j) {                   //refines the guesses
                    q = q * (z[j]- z[i]);
                }
            }
            deltaz[j] = -(fofx/q);              //set delta z to -r/q
            if(cabs(deltaz[j]) > cabs(zmax)) {  //if the new delta z is greater
                zmax = cabs(deltaz[j]);         //then old zmax set zmax to new delta z
            }
        }
        printf("iter %i\n", k+1);               //print out rest of guesses
        for (int j = 0; j < n; j++){
            z[j] = z[j] + deltaz[j];
            printf("z[%i] = %0.10f + %0.10f i\n", j, creal(z[j]), cimag(z[j]));
        }
        
        if(cabs(zmax) <= FLT_EPSILON) {         //stop iterating if zmax is less then
           break;                               //machine epsilon because nothing else can be done
        }
    }
    return 0;
}








