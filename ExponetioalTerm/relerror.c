//  Created by Emanuel Brici on 8/31/15.
//  Copyright (c) 2015 Emanuel Brici. All rights reserved.
//

#include <stdio.h>
#include <math.h>

//Function that computes factorial
double fac(int x) {
    if (x == 0 || x == 1){
        return 1;
    }
    return (x * fac(x - 1));
}

//Function that computes the upper bound error for n
double relerror(int n){
    for(int i = 1; i <= n ; i++){
        double x = 2.0/fac(i+1);
        double y = pow((log(2)/2),(i+1));
        double err = x * y;
        printf("n= %d------> error = %0.15e\n", i, err);
    }
    return 0;
}

int main(int argc, const char * argv[]) {
    relerror(15);
    return 0;
}