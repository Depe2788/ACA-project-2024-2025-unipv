#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define M 3
#define N 4
#define T 2
#define S 3

struct matrix {
        int nrows;
        int ncols;
        double *mat;
};
struct squareMatrix {
        int n;
        double *mat;
};
struct vector {
        int length; 
        double *vect;
};

void readMatrix(FILE *f, struct squareMatrix *matrix);
void initializeArray(struct vector *v);
void printArray(struct vector *v);
void initializeMatrix(double *mat, int nrows, int ncols);
void printMatrix(double *mat, int nrows, int ncols);
void matrixMul(struct matrix *matrix1, struct matrix *matrix2, struct matrix *matrix3);
void forwardSubstitution(struct squareMatrix *A, struct vector *b, struct vector *x);
void backwardSubstitution(struct squareMatrix *A, struct vector *b, struct vector *x);
void matrixInverse(struct squareMatrix *A, struct squareMatrix *inverse);
void matrixInversePivoting(struct squareMatrix *A, struct squareMatrix *inverse);


void readMatrix(FILE *f, struct squareMatrix *matrix){
        char buf[10]; 
        fgets(buf, sizeof(buf), f);
        sscanf(buf, "%i", &matrix->n);
        matrix->mat = (double *)malloc(matrix->n * sizeof(double));
        for(int i = 0; i < (matrix->n * matrix->n); i++){     
                fgets(buf, sizeof(buf),f);
                matrix->mat[i]=atof(buf);
        }
}

void initializeArray(struct vector *v)
{
        for (int i = 0; i < v->length; ++i)
                v->vect[i] = rand()%10;
}

void initializeMatrix(double *mat, int nrows, int ncols)
{
        for (int i = 0; i < nrows; i++){
                for(int j = 0; j < ncols; j++){
                        mat[ncols*i+j] = rand()%10;  //between 0 and 9
                } 
        }
}

void printArray(struct vector *v)
{
        for (int i = 0; i < v->length; ++i){
                printf("%0.1f\n", v->vect[i]);                
        }
}

void printMatrix(double *mat, int nrows, int ncols)
{
        for (int i = 0; i < nrows; i++){
                for(int j = 0; j < ncols; j++){
                        printf("%0.2f ", mat[ncols*i+j]);
                } 
                printf("\n");
        }
}

void matrixMul(struct matrix *matrix1, struct matrix *matrix2, struct matrix *matrix3){
        if (matrix1->ncols != matrix2->nrows) {
                printf("The matrix multiplication can't be done\n");
                free(matrix1->mat);
                free(matrix2->mat);
                exit(1);
        }
        matrix3->nrows = matrix1->nrows;
        matrix3->ncols = matrix2->ncols;
        matrix3->mat = (double *)malloc(matrix3->nrows * matrix3->ncols * sizeof(double));

        for (int i = 0; i < matrix1->nrows; i++){
                for (int j = 0; j < matrix2->ncols; j++){
                        matrix3->mat[matrix2->ncols*i+j] = 0;
                        for (int k = 0; k < matrix1->ncols; k++){
                                matrix3->mat[matrix2->ncols*i+j] += matrix1->mat[matrix1->ncols*i+k] * matrix2->mat[matrix2->ncols*k+j];
                        }
                }
        }
}

//linear systems Ax = b; A must be a square matrix non singular
//forward substitution method for lower triangular systems 
void forwardSubstitution(struct squareMatrix *A, struct vector *b, struct vector *x){
        for (int i = 0; i < A->n; i++){
                x->vect[i] = b->vect[i];
                for (int j = 0; j < i; j++){
                     x->vect[i] -= A->mat[(A->n) * i + j] * x->vect[j];
                }
                if (A->mat[(A->n) * i + i] == 0) {
                        printf("Errore: diagonal element 0 in the forward substitution (A singular).\n");
                        exit(1);
                }
                x->vect[i] /= A->mat[(A->n) * i + i];
        }
}

//linear systems Ax = b; A must be a square matrix non singular
//backward substitution method for upper triangular systems 
void backwardSubstitution(struct squareMatrix *A, struct vector *b, struct vector *x){
        for (int i = (A->n) - 1; i >= 0; i--){
                x->vect[i] = b->vect[i];
                for (int j = i+1; j < A->n; j++){
                    x->vect[i] -= A->mat[(A->n) * i + j] * x->vect[j];
                }
                if (A->mat[(A->n) * i + i] == 0) {
                        printf("Errore: diagonal element 0 in the backward substitution (A singular).\n");
                        exit(1);
                }
                x->vect[i] /= A->mat[(A->n) * i + i];
        }
}

void matrixInverse(struct squareMatrix *A, struct squareMatrix *inverse){
        inverse->n = A->n;
        inverse->mat = (double *)malloc(inverse->n * inverse->n * sizeof(double));

        //A=LU find L, U
        struct squareMatrix L, U;
        L.n = A->n;
        L.mat = (double *)malloc(L.n * L.n * sizeof(double));
        U.n = A->n; 
        U.mat = (double *)malloc(U.n * U.n * sizeof(double));

        //initialize L: identity matrix and U = A
        for (int i = 0; i < A->n; i++){
                for(int j = 0; j < A->n; j++){  
                        if(i == j){
                                L.mat[(L.n) * i + j] = 1;
                        } else {
                                L.mat[(L.n) * i + j] = 0;
                        } 
                        U.mat[(U.n) * i + j] = A->mat[(A->n) * i + j];
                }
        }

        for (int k = 0; k < (A->n) - 1; k++){
                // If the pivot is equal to zero LU factorization can't go on 
                if (U.mat[U.n * k + k] == 0) {
                        printf("Pivot equal to zero, LU factorization can't go on.\n");
                        free(L.mat);
                        free(U.mat);
                        return;  
                }
                for (int i = k+1; i < A->n; i++){
                        L.mat[L.n * i + k] = U.mat[U.n * i + k] / U.mat[U.n * k + k];
                        for (int s = k; s < A->n; s++){
                                U.mat[U.n * i + s] -= L.mat[L.n * i + k] * U.mat[U.n * k + s];
                        }
                }
        }
        /*
        //output: L, U   
        printf("Matrix L:\n");
        printMatrix(L.mat, L.n, L.n);
        printf("Matrix U:\n");
        printMatrix(U.mat, U.n, U.n); 
        */    

        struct vector y; 
        y.length = A->n;
        y.vect = (double *)malloc(y.length * sizeof(double));
        struct vector e;
        e.length = A->n;
        e.vect = (double *)malloc(e.length * sizeof(double));
        struct vector c;
        c.length = A->n;
        c.vect = (double *)malloc(c.length * sizeof(double));

        for (int i = 0; i < A->n; i++){
                for(int k = 0; k < e.length; k++){
                        e.vect[k]=0;
                }
                e.vect[i] = 1;
                forwardSubstitution(&L, &e, &y);
                backwardSubstitution(&U, &y, &c);
                //copy c in inverse
                for(int k = 0; k < A->n; k++){
                        inverse->mat[(A->n) * k + i] = c.vect[k];
                }
        }

        free(L.mat);
        free(U.mat);
        free(y.vect);
        free(e.vect);
        free(c.vect);
}

//compute the inverse with LU pivoting
void matrixInversePivoting(struct squareMatrix *A, struct squareMatrix *inverse){
        inverse->n = A->n;
        inverse->mat = (double *)malloc(inverse->n * inverse->n * sizeof(double));

        //PA=LU find L, U, P
        struct squareMatrix L, U, P;
        L.n = A->n;
        L.mat = (double *)malloc(L.n * L.n * sizeof(double));
        U.n = A->n; 
        U.mat = (double *)malloc(U.n * U.n * sizeof(double));
        P.n = A->n; 
        P.mat = (double *)malloc(P.n * P.n * sizeof(double));

        //initialize L and P (identity matrices) and U = A
        for (int i = 0; i < A->n; i++){
                for(int j = 0; j < A->n; j++){  
                        if(i == j){
                                L.mat[(L.n) * i + j] = 1;
                                P.mat[(P.n) * i + j] = 1;
                        } else {
                                L.mat[(L.n) * i + j] = 0;
                                P.mat[(P.n) * i + j] = 0;
                        } 
                        U.mat[(U.n) * i + j] = A->mat[(A->n) * i + j];
                }
        }

        double tmp;
        int maxIndex; 
        for (int k = 0; k < (A->n) - 1; k++){
                tmp = U.mat[(U.n)*k+k];
                maxIndex = k;
                for(int j = k; j < (A->n); j++){
                        if(fabs(U.mat[(U.n) * j + k]) > fabs(tmp)){
                                tmp = U.mat[(U.n)*j+k];
                                maxIndex = j; 
                        }
                }

                // Controllo pivot nullo
                if (fabs(tmp) < 1e-10) {
                        printf("Error: Pivot 0 or too small.\n");
                        exit(EXIT_FAILURE);
                }

                for(int s = k; s < A->n; s++){
                        tmp = U.mat[(U.n) * k + s];
                        U.mat[(U.n) * k + s] = U.mat[(U.n) * maxIndex + s]; 
                        U.mat[(U.n) * maxIndex + s] = tmp;
                }

                for(int s = 0; s < A->n; s++){
                        tmp = P.mat[(P.n) * k + s];
                        P.mat[(P.n) * k + s] = P.mat[(P.n) * maxIndex + s]; 
                        P.mat[(P.n) * maxIndex + s] = tmp;
                }

                if(k >= 1){
                        for(int s = 0; s < k; s++){
                                tmp = L.mat[(L.n) * k + s];
                                L.mat[(L.n) * k + s] = L.mat[(L.n) * maxIndex + s];
                                L.mat[(L.n) * maxIndex + s] = tmp;
                        }
                        
                }

                for (int i = k+1; i < A->n; i++){
                        L.mat[L.n * i + k] = U.mat[U.n * i + k] / U.mat[U.n * k + k];
                        for (int s = k; s < A->n; s++){
                                U.mat[U.n * i + s] = U.mat[U.n * i + s] - L.mat[L.n * i + k] * U.mat[U.n * k + s];
                        }
                }
        }
        
        //output: L, U, P
        printf("Matrix L:\n");
        printMatrix(L.mat, L.n, L.n);
        printf("Matrix U:\n");
        printMatrix(U.mat, U.n, U.n); 
        printf("Matrix P:\n");
        printMatrix(P.mat, P.n, P.n); 

        
        struct vector y; 
        y.length = A->n;
        y.vect = (double *)malloc(y.length * sizeof(double));
        struct vector pe;
        pe.length = A->n;
        pe.vect = (double *)malloc(pe.length * sizeof(double));
        struct vector c;
        c.length = A->n;
        c.vect = (double *)malloc(c.length * sizeof(double));

        for (int i = 0; i < A->n; i++){
                for(int k = 0; k < pe.length; k++){
                        pe.vect[k] = P.mat[P.n * k + i];
                }
                forwardSubstitution(&L, &pe, &y);
                backwardSubstitution(&U, &y, &c);
                //copy c in inverse
                for(int k = 0; k < A->n; k++){
                        inverse->mat[(A->n) * k + i] = c.vect[k];
                }
        }
  

        free(L.mat);
        free(U.mat);
        free(P.mat);
        free(y.vect);
        free(pe.vect);
        free(c.vect);
}

int main(int argc, char* argv[])
{
        clock_t timer; 

        struct matrix matrix1, matrix2;

        srand(time(NULL));
        //matrix1 initialization M x N
        matrix1.nrows = M; 
        matrix1.ncols = N;
        matrix1.mat = (double *)malloc(matrix1.nrows * matrix1.ncols * sizeof(double));
        initializeMatrix(matrix1.mat, matrix1.nrows, matrix1.ncols);
        puts("matrix1:");
        printMatrix(matrix1.mat, matrix1.nrows, matrix1.ncols);
        //matrix2 initialization N x P
        matrix2.nrows = N; 
        matrix2.ncols = T;
        matrix2.mat = (double *)malloc(matrix2.nrows * matrix2.ncols * sizeof(double));
        initializeMatrix(matrix2.mat, matrix2.nrows, matrix2.ncols);
        puts("matrix2:");
        printMatrix(matrix2.mat, matrix2.nrows, matrix2.ncols);

        //matrix multiplication result M x P
        struct matrix matrix3;
        timer = clock();
        matrixMul(&matrix1, &matrix2, &matrix3);
        timer = clock() - timer; 
        puts("matrix3:");
        printMatrix(matrix3.mat, matrix3.nrows, matrix3.ncols);
        printf("Time to compute matrix multiplication: %0.6f seconds\n", ((double)timer)/CLOCKS_PER_SEC);

        //inverse of a square matrix matrix4 read from file matrix4.txt
        struct squareMatrix matrix4, inverse;
        FILE *f; 
        f = fopen("matrix4.txt", "r");
        readMatrix(f, &matrix4);
        puts("matrix4:");
        printMatrix(matrix4.mat, matrix4.n, matrix4.n);
        timer = clock();
        matrixInversePivoting(&matrix4, &inverse);
        timer = clock() - timer; 
        puts("inverse:");
        printMatrix(inverse.mat, inverse.n, inverse.n);
        printf("Time to compute the inverse: %0.6f seconds\n", ((double)timer)/CLOCKS_PER_SEC);

        
        free(matrix1.mat);
        free(matrix2.mat);
        free(matrix3.mat);
        free(matrix4.mat);
        free(inverse.mat);

        return 0; 
}





        


