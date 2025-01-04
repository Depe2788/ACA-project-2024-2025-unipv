#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#define M 3
#define N 4
#define T 2

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
        matrix->mat = (double *)malloc(matrix->n * matrix->n * sizeof(double));
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

void transposeInPlace(struct squareMatrix *matrix) {
    double temp;
    for (int i = 0; i < matrix->n; i++) {
        for (int j = i + 1; j < matrix->n; j++) {
            temp = matrix->mat[matrix->n * i + j];
            matrix->mat[matrix->n * i + j] = matrix->mat[matrix->n * j + i];
            matrix->mat[matrix->n * j + i] = temp;
        }
    }
}

int main(int argc, char* argv[])
{
    //inverse of a square matrix matrix4 read from file matrix4.txt
    struct squareMatrix matrix4, inverse;
    struct matrix inversePart;
    
    FILE *f; 
    double timer;

    MPI_Init(&argc, &argv);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(my_rank == 0){
        f = fopen("matrix4.txt", "r");
        readMatrix(f, &matrix4);
        //printMatrix(matrix4.mat, matrix4.n, matrix4.n);
        inverse.n = matrix4.n;
        inverse.mat = (double *)malloc(inverse.n * inverse.n * sizeof(double));
        //printf("matrix4.n: %d\n", matrix4.n);
    }

    //share number of the elements of the inverse
    MPI_Bcast(&matrix4.n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    inversePart.ncols = matrix4.n;

    int *recvcounts;
    int *displs;

    if (my_rank == 0) {
            recvcounts = (int *)malloc(size * sizeof(int));
            displs = (int *)malloc(size * sizeof(int));

            int remaining = matrix4.n % size;
            for (int i = 0; i < size; i++) {
                    recvcounts[i] = (matrix4.n / size) * matrix4.n;
                    if (i < remaining) {
                            recvcounts[i] += matrix4.n;  
                    }
            }

            displs[0] = 0;
            for (int i = 1; i < size; i++) {
                    displs[i] = displs[i - 1] + recvcounts[i - 1];
            }
    }

    MPI_Scatter(recvcounts, 1, MPI_INT, &inversePart.nrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    inversePart.nrows /= inversePart.ncols; 
    printf("process %d, matrix1Part.nrows = %d\n", my_rank, inversePart.nrows);
    inversePart.mat = (double *)malloc(inversePart.nrows * inversePart.ncols * sizeof(double));
    int offset;
    MPI_Scatter(displs, 1, MPI_INT, &offset, 1, MPI_INT, 0, MPI_COMM_WORLD);
    offset /= inversePart.ncols;
    printf("process %d, offset = %d\n", my_rank, offset);
    
    struct squareMatrix L, U, P; 
    L.n = matrix4.n;
    U.n = matrix4.n; 
    P.n = matrix4.n; 
    L.mat = (double *)malloc(L.n * L.n * sizeof(double));
    U.mat = (double *)malloc(U.n * U.n * sizeof(double));
    P.mat = (double *)malloc(P.n * P.n * sizeof(double));

    struct vector y; 
    y.length = matrix4.n;
    y.vect = (double *)malloc(y.length * sizeof(double));
    struct vector pe;
    pe.length = matrix4.n;
    pe.vect = (double *)malloc(pe.length * sizeof(double));
    struct vector c;
    c.length = matrix4.n;
    c.vect = (double *)malloc(c.length * sizeof(double));

    if(my_rank == 0){
        timer = MPI_Wtime();
        //initialize L and P (identity matrices) and U = A
        for (int i = 0; i < matrix4.n; i++){
            for(int j = 0; j < matrix4.n; j++){  
                if(i == j){
                        L.mat[(L.n) * i + j] = 1;
                        P.mat[(P.n) * i + j] = 1;
                } else {
                        L.mat[(L.n) * i + j] = 0;
                        P.mat[(P.n) * i + j] = 0;
                } 
                U.mat[(U.n) * i + j] = matrix4.mat[(matrix4.n) * i + j];
            }
        }

        double tmp;
        int maxIndex; 
        for (int k = 0; k < (matrix4.n) - 1; k++){
            tmp = U.mat[(U.n)*k+k];
            maxIndex = k;
            for(int j = k; j < (matrix4.n); j++){
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

            for(int s = k; s < matrix4.n; s++){
                    tmp = U.mat[(U.n) * k + s];
                    U.mat[(U.n) * k + s] = U.mat[(U.n) * maxIndex + s]; 
                    U.mat[(U.n) * maxIndex + s] = tmp;
            }

            for(int s = 0; s < matrix4.n; s++){
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

            for (int i = k+1; i < matrix4.n; i++){
                    L.mat[L.n * i + k] = U.mat[U.n * i + k] / U.mat[U.n * k + k];
                    for (int s = k; s < matrix4.n; s++){
                            U.mat[U.n * i + s] = U.mat[U.n * i + s] - L.mat[L.n * i + k] * U.mat[U.n * k + s];
                    }
            }
        }
    }

    //forward L and U to other processes
    MPI_Bcast(L.mat, L.n * L.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(U.mat, U.n * U.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(P.mat, P.n * P.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < inversePart.nrows; i++){
        for(int k = 0; k < pe.length; k++){
            pe.vect[k] = P.mat[P.n * k + i + offset];
        }
        forwardSubstitution(&L, &pe, &y);
        backwardSubstitution(&U, &y, &c);
        //copy c in inversePart
        for(int k = 0; k < inversePart.ncols; k++){
            inversePart.mat[inversePart.ncols * i + k] = c.vect[k];         
        }
    }

    MPI_Gatherv(inversePart.mat, inversePart.nrows * inversePart.ncols, MPI_DOUBLE, inverse.mat, 
            recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(my_rank == 0){
        transposeInPlace(&inverse);
        timer = MPI_Wtime() - timer;
        printMatrix(inverse.mat, inverse.n, inverse.n);
        printf("Time to compute the inverse: %.6f\n", timer);
        free(inverse.mat);
    }

    free(L.mat);
    free(U.mat);
    free(P.mat);
    free(y.vect);
    free(pe.vect);
    free(c.vect);
    free(inversePart.mat);

    MPI_Finalize();

    return EXIT_SUCCESS;
}





        


