#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include<sys/time.h>

double get_time(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);

  return (tv.tv_sec) + 1.0e-6 * tv.tv_usec;
}

//this is the polynomial that is integrated. THe degree is specified by the user
//and (n+1)/2 is the most optimized degree of the integral
double function(double ex, int n){
    return pow(ex,n);
} //return exp(ex);}

//this is the Legendre polynomial and the degree is specified by the user
double p(int n, double x){
    if(n == 0.){//this is the case for n=0
        return 1.;
    } else if(n == 1.){
        return x;
    }else //this is the case for n>1
    return ((2.*(n-1.)+1)*x*p(n-1.,x)-(n-1.)*p(n-2.,x))/(n);
    //normally this recursive formula looks a bit different, but this one calculates n using n-1 and //n-2
}

int main(int argc, char **argv)
{
    int n;//this is the degree of the Legendre polynomial
    n = atoi(argv[1]);
    double a;//lower bound on the integral
    double b;//upper bound
    a = atof(argv[2]);
    b = atof(argv[3]);

    double *w;//these are the weights of the Legendre polynomials
    double *pPrime;//this is the derivative of p

    if(n>50 || n<=0){
        printf("That's not within the degree of accuracy\n");
        return 1;
    }

    // initializing MPI
    MPI_Init(&argc, &argv);
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const double t_start = get_time();

    int start = rank * (1000000 / size); // calculate the start and end values of v for this process
    int end = (rank + 1) * (1000000 / size);

    double *halfroot;
    double *halfrootp;
    int num_halfroots = 0; // variable to hold the total number of halfroots found by all ranks

    if (rank == 0){
        if(n % 2 == 1){//when n is odd, zero is always a root and it's never a root when n is even
            halfroot = (double*)malloc(((n+1)/2) * sizeof(double));
        }else{ 
            halfroot = (double*)malloc((n/2) * sizeof(double));
        }
    }

    if(n % 2 == 1){//when n is odd, zero is always a root and it's never a root when n is even
        halfrootp = (double*)malloc(((n+1)/2) * sizeof(double));
    }else{ 
        halfrootp = (double*)malloc((n/2) * sizeof(double));
    }

    w = (double*)malloc(n * sizeof(double)); 
    pPrime = (double*)malloc(n * sizeof(double));

    int i = -1;
    //below is the root calculator
    double root[n];//this stores all the roots
    for(int v = start; v<end; v++) {
        double x = v/1000000.0;
            if(p(n,x)*p(n,x+0.000001) < 0) {
                i++;
                halfrootp[i] = x;
        }
    }

    num_halfroots = i + 1; // store the number of halfroots found by this rank

    // Gather the number of halfroots found by each rank
    int *num_halfroots_list = (int*)malloc(size * sizeof(int));
    MPI_Allgather(&num_halfroots, 1, MPI_INT, num_halfroots_list, 1, MPI_INT, MPI_COMM_WORLD);

    // Calculate the displacement and recvcounts for the gather operation
    int *displs = (int*)malloc(size * sizeof(int));
    int *recvcounts = (int*)malloc(size * sizeof(int));
    int total_halfroots = 0;
    for(int j = 0; j < size; j++) {
        displs[j] = total_halfroots;
        recvcounts[j] = num_halfroots_list[j];
        total_halfroots += num_halfroots_list[j];
    }

    // Gather all the halfroots onto rank 0
    MPI_Gatherv(halfrootp, num_halfroots, MPI_DOUBLE, halfroot, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(rank == 0) { // print the halfroots only on rank 0
        /*for(int j = 0; j < total_halfroots; j++) {
            printf("half = %lf\n",halfroot[j]);
        }*/

        //this fills in the negative roots of the polynomial and 0 if n is odd
        if(n % 2 == 1){    
            for(int j=0; j<(n+1)/2; j++){//not filling in the last elements allows zero to be a root in //the odd array
                root[2*j+1] = -1*halfroot[j];
                root[2*j] = halfroot[j];
            }
        }
        else{//this is the even case
            for(int j=0; j<(n)/2; j++){
                root[2*j]=halfroot[j];
                root[2*j+1]=-1*halfroot[j];
                }
            }
    
    
        /*for(int j=0; j<n; j++){
            printf("root = %lf\n", root[j]);
        }*/

            //this calculates the derivative
        for(int l=0; l<=n-1;l++){//this is the derivative of the Legendre polynomial
        pPrime[l] = (root[l]*p(n,root[l])-p(n-1,root[l]))*(n/(root[l]*root[l]-1));
        //printf("%lf\n",pPrime[l]);
        }

        //this calculates the weights
        for(int h=0; h<=n-1;h++){
        w[h]=2/((1-(root[h]*root[h]))*(pPrime[h]*pPrime[h]));
        //printf("weight = %lf\n",w[h]);
        }

        double sum = 0.0;
        for(int i=0;i<=n-1;i++){//Legendre summation
           sum += w[i]*function((b-a)/2.*root[i]+(a+b)/2.,2.*n-1.);
        }
    
        sum = (b-a)/2.*sum;
        /*for(int g=0;g<=n-1;g++){
            printf("root = %lf\n",root[g]);
            printf("weight = %lf\n",w[g]);
        }*/
        const double t_end = get_time();
        printf("%lf\t%.2f\n", sum, t_end-t_start);

        free(halfroot);
        free(halfrootp);
        free(w);
        free(pPrime);

        }

    MPI_Finalize(); // Finalize MPI
    return 0;
}


