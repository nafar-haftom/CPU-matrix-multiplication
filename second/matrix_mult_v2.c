#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void fill(double* x, int n) {

    int i;
    for (i=0, n=n*n; i<n; i++, x++)
        *x = ((double) (1 + rand() % 12345)) / ((double) (1 + rand() % 6789));
}

void matrix_mult_index (int n, double* a, double* b, double* c) {
  int i, j, k;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      c[i*n+j] =0;
      for(k = 0; k < n; k++)
        c[i*n+j] += a[i*n+k] * b[k*n+j];
    }
}

void matrix_mult_ptr_reg (int n, double* a, double* b, double* c) {
    register double cij;
    register double *at, *bt;
    register int i, j, k;
    for (i=0; i<n; i++, a+=n)
        for (j = 0; j < n; j++, c++) {
            cij = 0;
            for(k = 0, at = a, bt = &b[j]; k < n; k++, at++, bt+=n)
                cij += *at * *bt;
            *c = cij;
        }
}

void matrix_mult_ptr_no_reg (int n, double* a, double* b, double* c) {
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
void transpose (int n, double* bt, double* b){
    register int i, j;
    register double temp;
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            temp = *(b + (i * n) + j);
            *(bt + (j * n) + i) = temp;
        }
    }
}

void matrix_mult_transpose(int n, double* a, double* b, double* c){
    double *B2  = (double*)_aligned_malloc(n * n * sizeof(double), 64 /*sizeof(double)*/);
    register double cij;
    register double *at, *bt,*bt1;
    register int i, j, k;
    transpose(n, B2,b);
    for (i=0 ,bt1=B2; i<n; i++, a+=n ,bt1=B2){
        for (j = 0; j < n; j++, c++ , bt1+=n) {
            cij = 0;
            for(k = 0, at = a, bt = bt1; k < n; k++, at++, bt++)
                cij += *at * *bt;
            *c = cij;
}   }  _aligned_free(B2);
}
void matrix_mult_block(int n, int block_size, double* a, double* b, double* c){
    double *B2  = (double*)_aligned_malloc(n * n * sizeof(double), 64 /*sizeof(double)*/);
    transpose(n, B2,b);
    register double cij;    
    register double *at, *bt, *at1, *bt1, *bt2, *at2 ,*bt3,*ct,*ct1,*ct2;
    register int i, j, k ,imax, jmax, kmax , i1 , j1 , k1;
    register int block_change= n * block_size;
    ct=c;
    for(i=0; i<n; i++) {
        for(j=0; j<n; j++, c++) {
            *c = 0;
}   }
    c=ct;
    ct1=c;
    ct = c;
    for (i=0 ,bt3=B2 ,ct2=c; i<n; i+=block_size, a+=block_change ,ct2+=block_change ,bt3=B2 ){

        imax = i + block_size > n ? n : i + block_size;

        for (j = 0 ,ct1=ct2; j < n; j+=block_size , bt3+=block_change , ct1+=block_size) {

            jmax = j + block_size > n ? n : j + block_size;

            for (k = 0, at2 = a, bt2 = bt3  ; k < n ; k+=block_size , at2+=block_size , bt2+=block_size ){

                kmax = k + block_size > n ? n : k + block_size;

                for (i1=i  ,bt1=bt2 , ct=ct1 ,at1=at2; i1<imax ; ++i1, at1+=n ,bt1=bt2 ,ct+=n  ){

                    for (j1 = j , c=ct ; j1 < jmax ; ++j1, c++ , bt1+=n) {
                        cij = 0;         
                        for (k1 = k, at = at1, bt = bt1; k1 < kmax ; ++k1, at++, bt++){                            
                            cij += *at * *bt;
                        }
                        *c += cij;
                    }
                }
            }
        }
    }
    _aligned_free(B2);
}
int main()
{
    //#pragma optimize("", on);
    clock_t t0, t1;
    int n, ref;

    do{
        printf("Input size of matrix, n = ");
        scanf("%d", &n);

        ref = 0;

        double *A  = (double*)_aligned_malloc(n * n * sizeof(double), 64 /*sizeof(double)*/); //  64 is cache line size
        double *B  = (double*)_aligned_malloc(n * n * sizeof(double), 64 /*sizeof(double)*/);
        double *C1 = (double*)_aligned_malloc(n * n * sizeof(double), 64 /*sizeof(double)*/);
        double *C2 = (double*)_aligned_malloc(n * n * sizeof(double), 64 /*sizeof(double)*/);

        if(A == NULL || B == NULL || C1 == NULL || C2 == NULL){
            printf("Memory Allocation Error\n\n");
            return(-1);
        }

        unsigned int seed = time(NULL);
        printf("\nseed = %u\n", seed);

        srand(seed);
        fill(A, n);
        fill(B, n);

        fflush(stdin);
        printf("\n\nDo you want to run matrix_mult_index (y/n)? ");
        if(getchar() == 'y'){
            ref = 1;
            t0 = clock();
            matrix_mult_index(n, A, B, C1);
            t1 = clock();
            printf("\n\t\t\tExecution time of matrix_mult_index = %0.2f s", (float)(t1-t0)/CLOCKS_PER_SEC);
        }

        fflush(stdin);
        printf("\n\nDo you want to run matrix_mult_ptr_reg (y/n)? ");
        if(getchar() == 'y'){
            ref++;
            t0 = clock();
            matrix_mult_ptr_reg(n, A, B, ref == 1 ? C1 : C2);
            t1 = clock();
            printf("\n\t\t\tExecution time of matrix_mult_ptr_reg = %0.2f s\n", (float)(t1-t0)/CLOCKS_PER_SEC);
        }

        fflush(stdin);
        printf("\n\nDo you want to run matrix_mult_ptr_no_reg (y/n)? ");
        if(getchar() == 'y'){
            if(++ref > 2) ref = 2;
            t0 = clock();
            matrix_mult_ptr_no_reg(n, A, B, ref == 1 ? C1 : C2);
            t1 = clock();
            printf("\n\t\t\tExecution time of matrix_mult_ptr_no_reg = %0.2f s", (float)(t1-t0)/CLOCKS_PER_SEC);
        }

        fflush(stdin);
        printf("\n\nDo you want to run matrix_mult_transpose (y/n)? ");
        if(getchar() == 'y'){
            if(++ref > 2) ref = 2;
            t0 = clock();
            matrix_mult_transpose(n, A, B, ref == 1 ? C1 : C2);
            t1 = clock();
            printf("\n\t\t\tExecution time of matrix_mult_transpose = %0.2f s", (float)(t1-t0)/CLOCKS_PER_SEC);
        }


        fflush(stdin);
        printf("\n\nDo you want to run matrix_mult_block (y/n)? ");
        if(getchar() == 'y'){
            if(++ref > 2) ref = 2;

            int block_size;
            printf("\n\tInput size of block = ");
            scanf("%d", &block_size);

            t0 = clock();
            matrix_mult_block(n, block_size, A, B, ref == 1 ? C1 : C2);
            t1 = clock();
            printf("\n\t\t\tExecution time of matrix_mult_block = %0.2f s", (float)(t1-t0)/CLOCKS_PER_SEC);
        }



        printf("\n\n\nEnd Of Execution\n\n");

        if(ref == 2){
            int i;
            double *c1, *c2;
            printf("\n\nStart of Compare: ");
            for(i=0, c1=C1, c2=C2, n=n*n; i<n; i++, c1++, c2++){
//              if(*c1 != *c2)
                if(fabs((*c1 - *c2) / *c1) > 1E-10)
                    break;
                if(i % (n/20) == 0)
                    printf(".");
            }

            if(i != n)
                printf(" Ooops, Error Found @ %d: %0.3f vs %0.3f\n\n",i, *c1, *c2);
            else
                printf(" OK, OK, Matrixes are equivalent.\n\n");
        }
        else
            printf("\n\nNo Compare due to No Reference or No Data.\n\n");

        _aligned_free(A);
        _aligned_free(B);
        _aligned_free(C1);
        _aligned_free(C2);

        fflush(stdin);
        printf("\n\nDo you want to continue (y/n)? ");

    } while(getchar() == 'y');

    return 0;
}
