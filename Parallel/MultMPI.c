#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define M 100
#define N 700
#define T 1000

#define MAX 15

struct matrix {
        int nrows;
        int ncols;
        double *mat;
};

void printMatrixFile(FILE *f, struct matrix *matrix);
void printMatrix(double *mat, int nrows, int ncols);
void initializeMatrix(double *mat, int nrows, int ncols);
void matrixMul(struct matrix *matrix1, struct matrix *matrix2, struct matrix *matrix3);

void printMatrixFile(FILE *f, struct matrix *matrix){
        for(int i = 0; i < matrix->nrows; i++){   
                for(int j = 0; j < matrix->ncols; j++){   
                        fprintf(f, "%.2f ", matrix->mat[matrix->ncols * i + j]);
                }
                fprintf(f, "\n");
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

void initializeMatrix(double *mat, int nrows, int ncols)
{
        for (int i = 0; i < nrows; i++){
                for(int j = 0; j < ncols; j++){
                        mat[ncols*i+j] = MAX * ((double)rand() / RAND_MAX);  //between 0 and MAX
                } 
        }
}

void matrixMul(struct matrix *matrix1, struct matrix *matrix2, struct matrix *matrix3){
        if (matrix1->ncols != matrix2->nrows) {
                printf("The matrix multiplication can't be done\n");
                free(matrix1->mat);
                free(matrix2->mat);
                exit(1);
        }

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
        double mul_timer, tot_timer;
        double mul_timer_max, tot_timer_max;

        FILE *f;

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
        matrix3.nrows = matrix1.nrows;
        matrix3.ncols = matrix2.ncols;

        struct matrix matrix3Part;
        matrix3Part.ncols = matrix2.ncols;
        
        if (my_rank == 0){

                srand(time(NULL));
               
                //matrix1 initialization M x N
                matrix1.mat = (double *)malloc(matrix1.nrows * matrix1.ncols * sizeof(double));
                initializeMatrix(matrix1.mat, matrix1.nrows, matrix1.ncols);
                f = fopen("matrix1.txt", "w");
                printMatrixFile(f, &matrix1);
                fclose(f);

                //matrix2 initialization N x P
                initializeMatrix(matrix2.mat, matrix2.nrows, matrix2.ncols);
                f = fopen("matrix2.txt", "w");
                printMatrixFile(f, &matrix2);
                fclose(f);
                
                //allocation of matrix3 (result of the mult) 
                matrix3.mat = (double *)malloc(matrix3.nrows * matrix3.ncols * sizeof(double));
        }
        
        int *sendcounts;
        int *displs;

        if (my_rank == 0) {
                sendcounts = (int *)malloc(size * sizeof(int));
                displs = (int *)malloc(size * sizeof(int));

                int remaining = matrix1.nrows % size;
                for (int i = 0; i < size; i++) {
                        sendcounts[i] = (matrix1.nrows / size) * matrix1.ncols;
                        if (i < remaining) {
                                sendcounts[i] += matrix1.ncols;  
                        }
                }

                displs[0] = 0;
                for (int i = 1; i < size; i++) {
                        displs[i] = displs[i - 1] + sendcounts[i - 1];
                }
        }
        
        MPI_Scatter(sendcounts, 1, MPI_INT, &matrix1Part.nrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
        matrix1Part.nrows /= matrix1Part.ncols; 
        printf("process %d, matrix1Part.nrows = %d\n", my_rank, matrix1Part.nrows);
        
        if (matrix1Part.nrows == 0){
                printf("More process than rows");
                MPI_Abort(MPI_COMM_WORLD, 1);
                MPI_Finalize();
                return 0;
        }

        matrix1Part.mat = (double *)malloc(matrix1Part.nrows * matrix1Part.ncols * sizeof(double));
        matrix3Part.nrows = matrix1Part.nrows;
        matrix3Part.mat = (double *)malloc(matrix3Part.nrows * matrix3Part.ncols * sizeof(double));
   
        tot_timer = MPI_Wtime();
        //send a number of rows of matrix1 to each process
        MPI_Scatterv(matrix1.mat, sendcounts, displs, MPI_DOUBLE, 
                matrix1Part.mat, matrix1Part.nrows * matrix1Part.ncols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        //broadcast matrix2 to all processes
        MPI_Bcast(matrix2.mat, matrix2.nrows * matrix2.ncols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        mul_timer = MPI_Wtime();
        matrixMul(&matrix1Part, &matrix2, &matrix3Part);
        mul_timer = MPI_Wtime() - mul_timer;
        
        if (my_rank == 0) {
                for (int i = 0; i < size; i++) {
                        sendcounts[i] = sendcounts[i] / matrix1.ncols * matrix3.ncols;
                }

                displs[0] = 0;
                for (int i = 1; i < size; i++) {
                        displs[i] = displs[i - 1] + sendcounts[i - 1];
                }
        }
        
        MPI_Gatherv(matrix3Part.mat, matrix3Part.nrows * matrix3Part.ncols, MPI_DOUBLE, 
                matrix3.mat, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        tot_timer = MPI_Wtime() - tot_timer;

        //printf("Process %d Time to compute matrix multiplication: %0.6f seconds\n", my_rank, mul_timer);
        //printf("Process %d Time to compute matrix multiplication and to transfer data: %0.6f seconds\n", my_rank, tot_timer);

        MPI_Reduce(&mul_timer, &mul_timer_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&tot_timer, &tot_timer_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        if (my_rank == 0){
                printf("Max process Time to compute matrix multiplication: %0.6f seconds\n", mul_timer_max);
                printf("Max process Time to compute matrix multiplication and to transfer data: %0.6f seconds\n", tot_timer_max);
                f = fopen("matrix3.txt", "w");
                printMatrixFile(f, &matrix3);
                fclose(f);
                free(matrix1.mat);
                free(matrix3.mat);
                free(displs);
                free(sendcounts);
        }
        free(matrix1Part.mat);
        free(matrix2.mat);
        free(matrix3Part.mat);
        
        MPI_Finalize();
 
        return EXIT_SUCCESS;
}



