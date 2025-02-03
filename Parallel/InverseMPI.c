#include "functions.h"

//inverse of a square matrix matrix4 read from file matrix4.txt
int main(int argc, char* argv[])
{
        struct squareMatrix matrix4, inverse;
        
        FILE *f; 
        double timer;

        MPI_Init(&argc, &argv);

        int my_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        //terminate useless processes
        if(my_rank >= matrix4.n){
                MPI_Finalize();
                exit(0);
        }

        if(my_rank == 0) {
                if(size > matrix4.n){
                        size = matrix4.n;
                }
                f = fopen("../Input/matrix4.txt", "r");
                //f = fopen("../Output/inverse.txt", "r");
                readSquareMatrixFile(f, &matrix4);  //only process 0 has the whole matrix4 
                fclose(f);
                inverse.n = matrix4.n;
                inverse.mat = (double *)malloc(inverse.n * inverse.n * sizeof(double)); //only process 0 has the whole inverse
        }

        //share number of the elements of the matrix to invert
        MPI_Bcast(&matrix4.n, 1, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
        if (my_rank == 0){
                timer = MPI_Wtime();
        }

        matrixInversePivoting(&matrix4, &inverse, my_rank, size, &timer);

        MPI_Barrier(MPI_COMM_WORLD);
        if(my_rank == 0){
                timer = MPI_Wtime() - timer;
                f = fopen("../Output/inverse.txt", "w");
                //f = fopen("../Input/original.txt", "w");
                printSquareMatrixFile(f, &inverse);
                fclose(f);
                printf("Time to compute the inverse: %.6f\n", timer);
                free(inverse.mat);
                free(matrix4.mat);
        }

        MPI_Finalize();

        return EXIT_SUCCESS;
}





        


