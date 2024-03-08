#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void fill(double* x, int n) {

    int i;
    for (i=0, n=n*n; i<n; i++, x++)
        *x = ((double) (1 + rand() % 12345)) / ((double) (1 + rand() % 6789));
}

void matrix_mult_0 (int n, double* a, double* b, double* c) {
  int i, j, k;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      c[i*n+j] =0;
      for(k = 0; k < n; k++)
        c[i*n+j] += a[i*n+k] * b[k*n+j];
    }
}

void matrix_mult_1 (int n, double* a, double* b, double* c) {
    double cij;
    double *at, *bt;
    int i, j, k;
    for (i=0; i<n; i++, a+=n)
        for (j = 0; j < n; j++, c++) {
            cij = 0;
            for(k = 0, at = a, bt = &b[j]; k < n; k++, at++, bt+=n)
                cij += *at * *bt;
            *c = cij;
        }
}

int main()
{
    int n;

    double *A, *B, *C1, *C2;

    printf("Input size of matrix, n = ");
    scanf("%d", &n);

    A =  (double*)malloc(n*n * sizeof(double));
    B =  (double*)malloc(n*n * sizeof(double));
    C1 = (double*)malloc(n*n * sizeof(double));
    C2 = (double*)malloc(n*n * sizeof(double));

    if(A == NULL || B == NULL || C1 == NULL || C2 == NULL){
        printf("Memory Allocation Error\n\n");
        return(-1);
    }

    srand(time(NULL));
    fill(A, n);
    fill(B, n);

    clock_t t0 = clock();
    matrix_mult_0(n, A, B, C1);
    clock_t t1 = clock();
    printf("\nExecution time of matrix_mult_0 = %0.3f s \n", (float)(t1-t0)/CLOCKS_PER_SEC);


    t0 = clock();
    matrix_mult_1(n, A, B, C2);
    t1 = clock();
    printf("\nExecution time of matrix_mult_1 = %0.3f s \n\n", (float)(t1-t0)/CLOCKS_PER_SEC);

    printf("End Of Execution\n\n\n\nStart of Compare: ");

    int i;
    for(i=0, n=n*n; i<n; i++){
        if(C1[i] != C2[i])
            break;
        if(i % (n/20) == 0)
            printf(".");
    }

    if(i != n)
        printf(" Ooops, Error Found @ %d: %f vs %f\n\n",i, C1[i], C2[i]);
    else
        printf(" OK, OK, Matrixes are equivalent.\n\n");

    free(A);
    free(B);
    free(C1);
    free(C2);

    return 0;
}
