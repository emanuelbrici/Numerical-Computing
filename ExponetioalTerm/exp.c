//  Created by Emanuel Brici on 8/31/15.
//  Copyright (c) 2015 Emanuel Brici. All rights reserved.
//

#include <stdio.h>
#include <math.h>

//Exponetial fuction that computes e^x with a Taylor series
float myexp(float x) {
    float z = x/M_LN2;              //equations broken down from assignmet 
    int m = round(z);               //given in class
    float M2 = ldexpf(1, m);
    float w = z - m;
    float u = w*log(2);
    
    const float u2 = u * u;         //Taylor series computation begins 
    float r = 1.0/5040;             //use of fmaf to keep the computation
    r = fmaf(r,u,1.0/720);          //of e^x as optimized as possible
    r = r * u2;
    r = fmaf(u, 1.0/120, r);
    r = (1.0/24 + r) * u2;
    r = fmaf(u, 1.0/6, r);
    r = (1.0/2 + r) * u2;
    r = r + 1.0 + u;
    r = M2*r;
    return r;
}

//function computes 30 values between -20 and 20
float comp() {
    const float x0 = -20;                //Range of numbers to use in myexp function
    const float x1 = 20;                   
    const float n = 30;
    const float dx = (x1 - x0)/(n-1);
    float x = x0;
    for (int i = 0; i < n; i++, x += dx) {
        float MYEXP = myexp(x);             //my exp function
        float EXP = exp(x);               
        float rerr = fabs((myexp(x)-exp(x))/exp(x)); //differnce between myexp and exp
        printf("%+0.9le %0.9le %0.9le %0.15e\n", x ,MYEXP, EXP, rerr);
    }
    return 0;
}

int main(int argc, const char * argv[]) {
    comp();
    return 0;
}

