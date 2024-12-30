#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define M 3
#define N 4
#define P 2
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

void initializeArray(struct vector *v);
void printArray(struct vector *v);
void initializeMatrix(double *mat, int nrows, int ncols);
void printMatrix(double *mat, int nrows, int ncols);
void matrixMul(struct matrix *matrix1, struct matrix *matrix2, struct matrix *matrix3);
void forwardSubstitution(struct squareMatrix *A, struct vector *b, struct vector *x);
void backwardSubstitution(struct squareMatrix *A, struct vector *b, struct vector *x);
void matrixInverse(struct squareMatrix *A, struct squareMatrix *inverse);


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
        //output: L, U   
        printf("Matrix L:\n");
        printMatrix(L.mat, L.n, L.n);
        printf("Matrix U:\n");
        printMatrix(U.mat, U.n, U.n);     

        //Ly=e1 find y with forward sub method
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


int main(int argc, char* argv[])
{
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
        matrix2.ncols = P;
        matrix2.mat = (double *)malloc(matrix2.nrows * matrix2.ncols * sizeof(double));
        initializeMatrix(matrix2.mat, matrix2.nrows, matrix2.ncols);
        puts("matrix2:");
        printMatrix(matrix2.mat, matrix2.nrows, matrix2.ncols);

        //matrix multiplication result M x P
        struct matrix matrix3;
        matrixMul(&matrix1, &matrix2, &matrix3);
        puts("matrix3:");
        printMatrix(matrix3.mat, matrix3.nrows, matrix3.ncols);

        //inverse of a square matrix matrix4
        struct squareMatrix matrix4, inverse;
        matrix4.n = S;
        matrix4.mat = (double *)malloc(matrix4.n * matrix4.n * sizeof(double));
        initializeMatrix(matrix4.mat, matrix4.n, matrix4.n);
        puts("matrix4:");
        printMatrix(matrix4.mat, matrix4.n, matrix4.n);
        
        inverse.n = matrix4.n;
        inverse.mat = (double *)malloc(matrix4.n * matrix4.n * sizeof(double));
        matrixInverse(&matrix4, &inverse);
        puts("inverse:");
        printMatrix(inverse.mat, inverse.n, inverse.n);
        
        free(matrix1.mat);
        free(matrix2.mat);
        free(matrix3.mat);
        free(matrix4.mat);
        free(inverse.mat);

        return 0; 
}





        


