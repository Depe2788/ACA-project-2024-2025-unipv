#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define M 8
#define N 4
#define T 2

struct matrix {
    int nrows;
    int ncols;
    double* mat;
};
struct squareMatrix {
    int n;
    double* mat;
};
struct vector {
    int length;
    double* vect;
};

void initializeArray(struct vector* v);
void printArray(struct vector* v);
void initializeMatrix(double* mat, int nrows, int ncols);
void printMatrix(double* mat, int nrows, int ncols);
void matrixMul(struct matrix* matrix1, struct matrix* matrix2, struct matrix* matrix3);
void forwardSubstitution(struct squareMatrix* A, struct vector* b, struct vector* x);
void backwardSubstitution(struct squareMatrix* A, struct vector* b, struct vector* x);
void matrixInversePivoting(struct squareMatrix* A, struct squareMatrix* inverse, int rank, int size);
void readMatrix(FILE* f, struct squareMatrix* matrix);
void printMatrixTranspose(double* mat, int nrows, int ncols);

void printMatrixTranspose(double* mat, int nrows, int ncols)
{
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            printf("%0.2f ", mat[ncols * j + i]);
        }
        printf("\n");
    }
}


void initializeArray(struct vector* v)
{
    for (int i = 0; i < v->length; ++i)
        v->vect[i] = rand() % 10;
}

void initializeMatrix(double* mat, int nrows, int ncols)
{
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            mat[ncols * i + j] = rand() % 10;  //between 0 and 9
        }
    }
}

void printArray(struct vector* v)
{
    for (int i = 0; i < v->length; ++i) {
        printf("%0.1f\n", v->vect[i]);
    }
}

void printMatrix(double* mat, int nrows, int ncols)
{
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            printf("%0.2f ", mat[ncols * i + j]);
        }
        printf("\n");
    }
}



void matrixMul(struct matrix* matrix1, struct matrix* matrix2, struct matrix* matrix3) {
    if (matrix1->ncols != matrix2->nrows) {
        printf("The matrix multiplication can't be done\n");
        free(matrix1->mat);
        free(matrix2->mat);
        exit(1);
    }
    matrix3->nrows = matrix1->nrows;
    matrix3->ncols = matrix2->ncols;
    matrix3->mat = (double*)malloc(matrix3->nrows * matrix3->ncols * sizeof(double));

    for (int i = 0; i < matrix1->nrows; i++) {
        for (int j = 0; j < matrix2->ncols; j++) {
            matrix3->mat[matrix2->ncols * i + j] = 0;
            for (int k = 0; k < matrix1->ncols; k++) {
                matrix3->mat[matrix2->ncols * i + j] += matrix1->mat[matrix1->ncols * i + k] * matrix2->mat[matrix2->ncols * k + j];
            }
        }
    }
}

//linear systems Ax = b; A must be a square matrix non singular
//forward substitution method for lower triangular systems
void forwardSubstitution(struct squareMatrix* A, struct vector* b, struct vector* x) {
    for (int i = 0; i < A->n; i++) {
        x->vect[i] = b->vect[i];
        for (int j = 0; j < i; j++) {
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
void backwardSubstitution(struct squareMatrix* A, struct vector* b, struct vector* x) {
    for (int i = (A->n) - 1; i >= 0; i--) {
        x->vect[i] = b->vect[i];
        for (int j = i + 1; j < A->n; j++) {
            x->vect[i] -= A->mat[(A->n) * i + j] * x->vect[j];
        }
        if (A->mat[(A->n) * i + i] == 0) {
            printf("Errore: diagonal element 0 in the backward substitution (A singular).\n");
            exit(1);
        }
        x->vect[i] /= A->mat[(A->n) * i + i];
    }
}

//compute the inverse with LU pivoting
void matrixInversePivoting(struct squareMatrix* A, struct squareMatrix* inverse, int rank, int size) {
    inverse->n = A->n;
    inverse->mat = (double*)malloc(inverse->n * inverse->n * sizeof(double));

    //PA=LU find L, U, P (only rank 0)
    struct squareMatrix L, U, P;
    L.n = A->n;
    L.mat = (double*)malloc(L.n * L.n * sizeof(double));
    U.n = A->n;
    U.mat = (double*)malloc(U.n * U.n * sizeof(double));
    P.n = A->n;
    P.mat = (double*)malloc(P.n * P.n * sizeof(double));

    if (rank == 0) {
        //Intialization of L, P as identity matrix and U = A
        for (int i = 0; i < A->n; i++) {
            for (int j = 0; j < A->n; j++) {
                if (i == j) {
                    L.mat[(L.n) * i + j] = 1;
                    P.mat[(P.n) * i + j] = 1;
                }
                else {
                    L.mat[(L.n) * i + j] = 0;
                    P.mat[(P.n) * i + j] = 0;
                }
                U.mat[(U.n) * i + j] = A->mat[(A->n) * i + j];
            }
        }

        //LU decomposition with pivoting (solo rank 0)
        double tmp;
        int maxIndex;
        for (int k = 0; k < (A->n) - 1; k++) {
            tmp = U.mat[(U.n) * k + k];
            maxIndex = k;
            for (int j = k; j < (A->n); j++) {
                if (fabs(U.mat[(U.n) * j + k]) > fabs(tmp)) {
                    tmp = U.mat[(U.n) * j + k];
                    maxIndex = j;
                }
            }

            // check pivot = 0
            if (fabs(tmp) < 1e-10) {
                printf("Error: Pivot 0 or too small.\n");
                exit(EXIT_FAILURE);
            }

            // row exchange between U and P
            for (int s = k; s < A->n; s++) {
                tmp = U.mat[(U.n) * k + s];
                U.mat[(U.n) * k + s] = U.mat[(U.n) * maxIndex + s];
                U.mat[(U.n) * maxIndex + s] = tmp;
            }

            for (int s = 0; s < A->n; s++) {
                tmp = P.mat[(P.n) * k + s];
                P.mat[(P.n) * k + s] = P.mat[(P.n) * maxIndex + s];
                P.mat[(P.n) * maxIndex + s] = tmp;
            }

            if (k >= 1) {
                for (int s = 0; s < k; s++) {
                    tmp = L.mat[(L.n) * k + s];
                    L.mat[(L.n) * k + s] = L.mat[(L.n) * maxIndex + s];
                    L.mat[(L.n) * maxIndex + s] = tmp;
                }
            }

            // Refresh of L and U
            for (int i = k + 1; i < A->n; i++) {
                L.mat[L.n * i + k] = U.mat[U.n * i + k] / U.mat[U.n * k + k];
                for (int s = k; s < A->n; s++) {
                    U.mat[U.n * i + s] = U.mat[U.n * i + s] - L.mat[L.n * i + k] * U.mat[U.n * k + s];
                }
            }
        }
    }
     

    MPI_Bcast(L.mat, L.n * L.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(U.mat, U.n * U.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(P.mat, P.n * P.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //Starting of the parallel inversion
    struct vector y, pe, c;
    y.length = A->n;
    pe.length = A->n;
    c.length = A->n;
    y.vect = (double*)malloc(y.length * sizeof(double));
    pe.vect = (double*)malloc(pe.length * sizeof(double));
    c.vect = (double*)malloc(c.length * sizeof(double));

    // Split column between processes
    int colsPerProcess = A->n / size;
    int startCol = rank * colsPerProcess;
    int endCol = (rank + 1) * colsPerProcess;
   
    struct matrix matrixTemp;
    matrixTemp.ncols = colsPerProcess;  
    matrixTemp.nrows = A->n;
    matrixTemp.mat = (double*)malloc(matrixTemp.nrows * matrixTemp.ncols * sizeof(double));

    int index=0;
    // Calcolo dell'inversa per le colonne locali
    for (int i = 0; i < colsPerProcess; i++) {
      int globalCol = rank * colsPerProcess + i;
      for (int k = 0; k < pe.length; k++) {
          pe.vect[k] = P.mat[P.n * k + globalCol];
      }
     
      forwardSubstitution(&L, &pe, &y);
      backwardSubstitution(&U, &y, &c);
     
      for (int k = 0; k < A->n; k++) {
        matrixTemp.mat[index++] = c.vect[k]; // Usa l'indice per salvare i valori in modo sequenziale
        printf("VAL: %f\n", c.vect[k]);
      }
    }
   
    //Intermediate results collection
    MPI_Gather(matrixTemp.mat, matrixTemp.ncols * matrixTemp.nrows, MPI_DOUBLE, inverse->mat, matrixTemp.ncols * matrixTemp.nrows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Free allocated memory
    free(y.vect);
    free(pe.vect);
    free(c.vect);
   
    //free(localA);
    if (rank == 0) {
        free(L.mat);
        free(U.mat);
        free(P.mat);
    }
}


void readMatrix(FILE* f, struct squareMatrix* matrix) {
    char buf[10];
    fgets(buf, sizeof(buf), f);
    sscanf(buf, "%i", &matrix->n);
    matrix->mat = (double*)malloc(matrix->n *matrix->n* sizeof(double));
    for (int i = 0; i < (matrix->n * matrix->n); i++) {
        fgets(buf, sizeof(buf), f);
        matrix->mat[i] = atof(buf);
    }
}

int main(int argc, char* argv[])
{
    clock_t timer;

    MPI_Init(&argc, &argv);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    struct matrix matrix1;
    matrix1.nrows = M;
    matrix1.ncols = N;

    struct matrix matrix2;
    matrix2.nrows = N;
    matrix2.ncols = T;
    matrix2.mat = (double*)malloc(matrix2.nrows * matrix2.ncols * sizeof(double));

    struct matrix matrix1Part;
    matrix1Part.ncols = matrix1.ncols;

    struct matrix matrix3;
    matrix3.nrows = M;
    matrix3.ncols = T;

    struct matrix matrix3Part;
    matrix3Part.ncols = T;

    struct squareMatrix matrix4, inverse;
    FILE* f;

    if (my_rank == 0) {

        srand(time(NULL));

        //matrix1 initialization M x N
        matrix1.mat = (double*)malloc(matrix1.nrows * matrix1.ncols * sizeof(double));
        initializeMatrix(matrix1.mat, matrix1.nrows, matrix1.ncols);
        puts("matrix1:");
        printMatrix(matrix1.mat, matrix1.nrows, matrix1.ncols);

        //matrix2 initialization N x P
        initializeMatrix(matrix2.mat, matrix2.nrows, matrix2.ncols);
        puts("matrix2:");
        printMatrix(matrix2.mat, matrix2.nrows, matrix2.ncols);

        // Allcoation of Matrix4 for inversion
        f = fopen("matrix4.txt", "r");
        readMatrix(f, &matrix4);
        puts("matrix4:");
        printMatrix(matrix4.mat, matrix4.n, matrix4.n);

        if (matrix1.nrows % size || matrix4.n % size) {
            printf("rows can't be splitted among the processes\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        matrix1Part.nrows = matrix1.nrows / size;
        printf("process %d send %d\n", my_rank, matrix1Part.nrows);

        //allocation of matrix3 (result of the mult)
        matrix3.mat = (double*)malloc(matrix3.nrows * matrix3.ncols * sizeof(double));

    }
   
    MPI_Bcast(&matrix4.n, 1, MPI_INT, 0, MPI_COMM_WORLD);
   
    //Inversion----------------------------------------------------
    clock_t timer_inv;
    timer_inv = clock();
    matrixInversePivoting(&matrix4, &inverse, my_rank, size);
    timer_inv = clock() - timer_inv;

    //send the number of rows manage by each process so that they can allocate their matrix1Part
    MPI_Bcast(&matrix1Part.nrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    matrix1Part.mat = (double*)malloc(matrix1Part.nrows * matrix1Part.ncols * sizeof(double));
    matrix3Part.nrows = matrix1Part.nrows;
    matrix3Part.mat = (double*)malloc(matrix3Part.nrows * matrix3Part.ncols * sizeof(double));

    if (my_rank != 0) {
        printf("process %d receive %d\n", my_rank, matrix1Part.nrows);
    }

    //send a number of rows of matrix1 to each process
    MPI_Scatter(matrix1.mat, matrix1Part.nrows * matrix1Part.ncols, MPI_DOUBLE,
        matrix1Part.mat, matrix1Part.nrows * matrix1Part.ncols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    printf("process %d, my matrix1Part is: \n", my_rank);
    printMatrix(matrix1Part.mat, matrix1Part.nrows, matrix1Part.ncols);

    //broadcast matrix2 to all processes
    MPI_Bcast(matrix2.mat, matrix2.nrows * matrix2.ncols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    timer = clock();
    //multiplication
    matrixMul(&matrix1Part, &matrix2, &matrix3Part);
    timer = clock() - timer;

    MPI_Gather(matrix3Part.mat, matrix3Part.nrows * matrix3Part.ncols, MPI_DOUBLE,
        matrix3.mat, matrix3Part.nrows * matrix3Part.ncols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (my_rank == 0) {
        printf("process %d, my matrix3 is: \n", my_rank);
        printMatrix(matrix3.mat, matrix3.nrows, matrix3.ncols);
        printf("Time to compute matrix multiplication: %0.6f seconds\n", ((double)timer) / CLOCKS_PER_SEC);
        free(matrix1.mat);

        puts("inverse:");
        printMatrixTranspose(inverse.mat, inverse.n,inverse.n);
        //printMatrix(inverse.mat, inverse.n, inverse.n);
        printf("Time to compute the inverse: %0.6f seconds\n", ((double)timer_inv) / CLOCKS_PER_SEC);
        free(matrix4.mat);
        free(inverse.mat);
    }
    free(matrix1Part.mat);
    free(matrix2.mat);


    MPI_Finalize();

    return EXIT_SUCCESS;
}