#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define M 2000
#define N 1000
#define T 1000

#define MAX 15

struct matrix {
        int nrows;
        int ncols;
        double *mat;
};
struct squareMatrix {
    int n;
    double* mat;
};
struct vector {
    int length;
    double* vect;
};

//general
void initializeMatrix(double *mat, int nrows, int ncols);
void printMatrix(double *mat, int nrows, int ncols);
void printMatrixFile(FILE *f, double *mat, int nrows, int ncols);
void readMatrixFile(FILE *f, double **mat, int *nrows, int *ncols);
void readSquareMatrixFile(FILE *f, struct squareMatrix *matrix);
void printSquareMatrixFile(FILE *f, struct squareMatrix *matrix);
void transposeInPlace(struct squareMatrix *matrix);
void printIntegerVector(int *v, int n);
void printVectorStruct(struct vector *v);

//multiplication
void matrixMul(struct matrix *matrix1, struct matrix *matrix2, struct matrix *matrix3);

//inverse
void forwardSubstitution(struct squareMatrix *A, struct vector *b);
void backwardSubstitution(struct squareMatrix *A, struct vector *b);
void matrixInversePivoting(struct squareMatrix *A, struct squareMatrix *inverse, int my_rank, int size);
void matrixInversePivotingImproved(struct squareMatrix *A, struct squareMatrix *inverse, int my_rank, int size);



