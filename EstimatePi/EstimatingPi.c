//
//  EstimatingPi.c
//  EstimatingPi
//
//  Created by Emanuel Brici on 10/26/15.
//  Copyright © 2015 Emanuel Brici. All rights reserved.
//

#include <stdio.h>
#include <assert.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

//Function that will be used
long double f(long double x) {
    return 4.0L/(1.0L + x*x);
}

//Implementation of Trapezoid Rule
long double trap(long double a, long double b, int n) {
    assert(a < b);
    const long double h = (b - a)/n;  //find with that will be used
    long double sum = 0.0;
    long double x = a + h;
    for (int i = 1; i < n; i++, x += h) // loop for every iteration starting at f1
        sum += f(x);
    return h*((f(a) + f(b))/2 + sum);
}

//Implementation of Simpsons 1/3 rule
long double simpsons(long double a, long double b, int n) {
    assert(n % 2 == 0);
    assert(a < b);
    const long double h = (b - a)/n;  //find with that will be used
    long double sum1 = 0.0;
    for (int i = 1; i < n; i += 2) // loop for every second iteration starting at f1
        sum1 += f(a + i*h);
    long double sum2 = 0.0;
    for (int i = 2; i < n; i += 2) // loop for every second iteration starting at f2
        sum2 += f(a + i*h);
    return h/3*((f(a) + f(b)) + 4*sum1 + 2*sum2);
}

//Implementation of Simpsons 3/8 rule
long double simpsons38(long double a, long double b, int n) {
    assert(a < b);
    const long double h = (b-a)/n;  //find with that will be used
    long double sum1 = 0.0;
    for (int i = 1; i < n; i += 3) // loop for every third iteration starting at f1
        sum1 += f(a + i*h);
    long double sum2 = 0.0;
    for (int i = 2; i < n; i += 3) // loop for every third iteration starting at f2
        sum2 += f(a + i*h);
    long double sum3 = 0.0;
    for (int i = 3; i < n; i += 3) // loop for every third iteration starting at f3
        sum3 += f(a + i*h);
    return ((3*h)/8)*((f(a) + f(b)) + 3*sum1 + 3*sum2 + 2*sum3); //result of
}

//Implementation of Booles rule
long double booles(long double a, long double b, int n) {
    assert(a < b);
    const long double h = (b-a)/n; //find with that will be used
    long double sum1 = 0.0;
    for (int i = 4; i < n; i += 4) //loop for every forth iteration starting at f4
        sum1 += 14*(f(a + i*h));
    long double sum2 = 0.0;
    for (int i = 1; i < n; i += 2) //loop for every second iteration starting at f1
        sum2 += 32*(f(a + i*h));
    long double sum3 = 0.0;
    for (int i = 2; i < n; i += 4) //loop for every forth iteration statinf at f2
        sum3 += 12*(f(a + i*h));
    return ((2*h)/45)*((7*(f(a) + f(b))) + sum1 + sum2 + sum3);
}

int main(void) {
    const unsigned maxn = 800000;
    printf("|     N | Trapezoid Error  |   Simpsons1/3 Error|   Simpsons3/8 Error|   Boole Error     |\n");
    for (unsigned n = 12; n <= maxn; n *= 2) {
        const long double I_trap = trap(0, 1, n);
        const long double error_trap = fabsl(M_PI - I_trap);
        const long double I_simp = simpsons(0, 1, n);
        const long double error_simp = fabsl(M_PI - I_simp);
        const long double I_simp38 = simpsons38(0, 1, n);
        const long double error_simp38 = fabsl(M_PI - I_simp38);
        const long double I_booles = booles(0, 1, n);
        const long double error_booles = fabsl(M_PI - I_booles);
        printf("|%6d | %0.10Le |   %0.10Le |   %0.10Le |   %0.10Le|\n", n,
               error_trap, error_simp, error_simp38, error_booles);
    }
    return 0;
}
