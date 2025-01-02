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
        matrix2.mat = (double *)malloc(matrix2.nrows * matrix2.ncols * sizeof(double));

        struct matrix matrix1Part;
        matrix1Part.ncols = matrix1.ncols;

        struct matrix matrix3;
        matrix3.nrows = M;
        matrix3.ncols = T;

        struct matrix matrix3Part;
        matrix3Part.ncols = T;
        
        if (my_rank == 0){

                srand(time(NULL));
               
                //matrix1 initialization M x N
                matrix1.mat = (double *)malloc(matrix1.nrows * matrix1.ncols * sizeof(double));
                initializeMatrix(matrix1.mat, matrix1.nrows, matrix1.ncols);
                puts("matrix1:");
                printMatrix(matrix1.mat, matrix1.nrows, matrix1.ncols);

                //matrix2 initialization N x P
                initializeMatrix(matrix2.mat, matrix2.nrows, matrix2.ncols);
                puts("matrix2:");
                printMatrix(matrix2.mat, matrix2.nrows, matrix2.ncols);

                if(matrix1.nrows % size){
                        printf("rows can't be splitted among the processes\n");
                        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                } 

                matrix1Part.nrows = matrix1.nrows / size;
                printf("process %d send %d\n", my_rank, matrix1Part.nrows);

                //allocation of matrix3 (result of the mult) 
                matrix3.mat = (double *)malloc(matrix3.nrows * matrix3.ncols * sizeof(double));
        }

        //send the number of rows manage by each process so that they can allocate their matrix1Part
        MPI_Bcast(&matrix1Part.nrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
        matrix1Part.mat = (double *)malloc(matrix1Part.nrows * matrix1Part.ncols * sizeof(double));
        matrix3Part.nrows = matrix1Part.nrows;
        matrix3Part.mat = (double *)malloc(matrix3Part.nrows * matrix3Part.ncols * sizeof(double));

        if(my_rank != 0){
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
        
        if (my_rank == 0){
                printf("process %d, my matrix3 is: \n", my_rank);
                printMatrix(matrix3.mat, matrix3.nrows, matrix3.ncols);
                printf("Time to compute matrix multiplication: %0.6f seconds\n", ((double)timer)/CLOCKS_PER_SEC);
                free(matrix1.mat);
        }
        free(matrix1Part.mat);
        free(matrix2.mat);

        MPI_Finalize();
 
        return EXIT_SUCCESS;
}



