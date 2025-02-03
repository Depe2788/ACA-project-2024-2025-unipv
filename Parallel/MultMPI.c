#include "functions.h"

int main(int argc, char* argv[])
{
        double mul_timer, tot_timer;

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

        int *sendcounts;
        int *displs;

        //terminate useless processes
        if(my_rank >= matrix1.nrows){
                MPI_Finalize();
                exit(0);
        }
        
        if (my_rank == 0){

                srand(time(NULL));
               
                //matrix1 initialization M x N
                matrix1.mat = (double *)malloc(matrix1.nrows * matrix1.ncols * sizeof(double));
                initializeMatrix(matrix1.mat, matrix1.nrows, matrix1.ncols);
                f = fopen("../Output/matrix1.txt", "w");
                printMatrixFile(f, matrix1.mat, matrix1.nrows, matrix1.ncols);
                fclose(f);

                //matrix2 initialization N x P
                initializeMatrix(matrix2.mat, matrix2.nrows, matrix2.ncols);
                f = fopen("../Output/matrix2.txt", "w");
                printMatrixFile(f, matrix2.mat, matrix2.nrows, matrix2.ncols);
                fclose(f);
                
                //allocation of matrix3 (result of the mult) 
                matrix3.mat = (double *)malloc(matrix3.nrows * matrix3.ncols * sizeof(double));

                if(size > matrix1.nrows){
                        size = matrix1.nrows;
                }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        if(my_rank == 0){
                tot_timer = MPI_Wtime();

                //sendcounts is the vector containing the number of elements of matrix1Part of each process
                sendcounts = (int *)malloc(size * sizeof(int));
                displs = (int *)malloc(size * sizeof(int));

                int remaining = matrix1.nrows % size;
                int div = matrix1.nrows / size;
                for (int i = 0; i < size; i++) {
                        sendcounts[i] = div * matrix1.ncols;
                        if (i < remaining) {
                                sendcounts[i] += matrix1.ncols;  
                        }
                }

                displs[0] = 0;
                for (int i = 1; i < size; i++) {
                        displs[i] = displs[i - 1] + sendcounts[i - 1];
                }
        }

        //MPI_Scatter ensures that data is sent from 0 only once ready, other processes are so implicitly blocked   
        MPI_Scatter(sendcounts, 1, MPI_INT, &matrix1Part.nrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
        matrix1Part.nrows /= matrix1Part.ncols; 
        //printf("process %d, matrix1Part.nrows = %d\n", my_rank, matrix1Part.nrows);

        matrix1Part.mat = (double *)malloc(matrix1Part.nrows * matrix1Part.ncols * sizeof(double));
        matrix3Part.nrows = matrix1Part.nrows;
        matrix3Part.mat = (double *)malloc(matrix3Part.nrows * matrix3Part.ncols * sizeof(double));
   
        //send a variable number of rows of matrix1 to each process
        MPI_Scatterv(matrix1.mat, sendcounts, displs, MPI_DOUBLE, 
                matrix1Part.mat, matrix1Part.nrows * matrix1Part.ncols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        //broadcast matrix2 to all processes
        MPI_Bcast(matrix2.mat, matrix2.nrows * matrix2.ncols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        matrixMul(&matrix1Part, &matrix2, &matrix3Part);
        
        if (my_rank == 0) {
                for (int i = 0; i < size; i++) {
                        sendcounts[i] = sendcounts[i] / matrix1.ncols * matrix3.ncols;
                }

                displs[0] = 0;
                for (int i = 1; i < size; i++) {
                        displs[i] = displs[i - 1] + sendcounts[i - 1];
                }
        }
        
        //MPI_Barrier(MPI_COMM_WORLD); not necessary
        //only when all processes call MPI_Gatherv data are united, so function blocking 
        MPI_Gatherv(matrix3Part.mat, matrix3Part.nrows * matrix3Part.ncols, MPI_DOUBLE, 
                matrix3.mat, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
        if (my_rank == 0){
                tot_timer = MPI_Wtime() - tot_timer;
                printf("Time to compute matrix multiplication and to transfer data: %0.6f seconds\n", tot_timer);
                f = fopen("../Output/matrix3.txt", "w");
                printMatrixFile(f, matrix3.mat, matrix3.nrows, matrix3.ncols);
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



