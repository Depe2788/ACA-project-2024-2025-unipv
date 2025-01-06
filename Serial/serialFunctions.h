#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define M 3
#define N 4
#define T 2
#define MAX 9

#define FILE_ERROR 1

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

//general
void readMatrixFile(FILE *f, double **mat, int *nrows, int *ncols);
void initializeMatrix(double *mat, int nrows, int ncols);
void printMatrix(double *mat, int nrows, int ncols);
void printMatrixFile(FILE *f, double *mat, int nrows, int ncols);

//multiplication
void matrixMul(struct matrix *matrix1, struct matrix *matrix2, struct matrix *matrix3);

//inverse
void forwardSubstitution(struct squareMatrix *A, struct vector *b, struct vector *x);
void backwardSubstitution(struct squareMatrix *A, struct vector *b, struct vector *x);
void matrixInversePivoting(struct squareMatrix *A, struct squareMatrix *inverse);
