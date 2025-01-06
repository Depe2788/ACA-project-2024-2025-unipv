#include "functions.h"

//compute the matrix multiplication
int main(int argc, char* argv[])
{
        clock_t timer; 

        struct matrix matrix1, matrix2;
        FILE *f;

        srand(time(NULL));
        //matrix1 initialization M x N
        matrix1.nrows = M; 
        matrix1.ncols = N;
        matrix1.mat = (double *)malloc(matrix1.nrows * matrix1.ncols * sizeof(double));
        initializeMatrix(matrix1.mat, matrix1.nrows, matrix1.ncols);
        f = fopen("../Output/matrix1.txt", "w");
        printMatrixFile(f, matrix1.mat, matrix1.nrows, matrix1.ncols);
        fclose(f);
        //matrix2 initialization N x P
        matrix2.nrows = N; 
        matrix2.ncols = T;
        matrix2.mat = (double *)malloc(matrix2.nrows * matrix2.ncols * sizeof(double));
        initializeMatrix(matrix2.mat, matrix2.nrows, matrix2.ncols);
        f = fopen("../Output/matrix2.txt", "w");
        printMatrixFile(f, matrix2.mat, matrix2.nrows, matrix2.ncols);
        fclose(f);

        //matrix multiplication result M x P
        struct matrix matrix3;
        timer = clock();
        matrixMul(&matrix1, &matrix2, &matrix3);
        timer = clock() - timer; 
        f = fopen("../Output/matrix3.txt", "w");
        printMatrixFile(f, matrix3.mat, matrix3.nrows, matrix3.ncols);
        fclose(f);
        printf("Time to compute matrix multiplication: %0.6f seconds\n", ((double)timer)/CLOCKS_PER_SEC);
        
        free(matrix1.mat);
        free(matrix2.mat);
        free(matrix3.mat);

        return 0; 
}





        


