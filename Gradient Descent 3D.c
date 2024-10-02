#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//function to minimize
double function(double x, double y){
    return x*x*x-8*x*x+y*y-y;
}

//derivative limit = h
double h = 0.00000001;

//derivative in the x direction
double xderivative(double x, double y){
    return (function(x+h,y)-function(x-h,y))/(2*h);
}

double yderivative(double x, double y){
    return (function(x,y+h)-function(x,y-h))/(2*h);
}

int main()
{
    //number of iterations and step size
    double step = 0.01;
    int iter = 1000;

    double *xcoords;
    double *ycoords;

    //holds the x values
    xcoords = (double*)malloc((iter) * sizeof(double));
    xcoords[0] = 8.3;

    //holds the y values
    ycoords = (double*)malloc((iter) * sizeof(double));
    ycoords[0] = 0.3;

    //main loop
    for(int i=0; i<iter; i++){
        //updates the function values and goes in the opposite direction of the gradient
        function(xcoords[i],ycoords[i]);
        xcoords[i+1] = xcoords[i] - step*xderivative(xcoords[i], ycoords[i]);
        ycoords[i+1] = ycoords[i] - step*yderivative(xcoords[i], ycoords[i]);
    }

    printf("x = %lf\n y = %lf\n z = %lf", xcoords[iter], ycoords[iter], function(xcoords[iter],ycoords[iter]));

    free(xcoords);
    free(ycoords);

    return 0;
}
