#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define M 10
#define N 8
#define T 10
#define MAX 15

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
void initializeMatrix(double *mat, int nrows, int ncols);
void printMatrix(double *mat, int nrows, int ncols);
void printMatrixFile(FILE *f, double *mat, int nrows, int ncols);
void printSquareMatrixFile(FILE *f, struct squareMatrix *matrix);
void readSquareMatrixFile(FILE *f, struct squareMatrix *matrix);

//multiplication
void matrixMul(struct matrix *matrix1, struct matrix *matrix2, struct matrix *matrix3);

//inverse
void forwardSubstitution(struct squareMatrix *A, struct vector *b);
void backwardSubstitution(struct squareMatrix *A, struct vector *b);
void matrixInversePivoting(struct squareMatrix *A, struct squareMatrix *inverse);
