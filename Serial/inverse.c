#include "functions.h"

//compute the inverse of a square matrix matrix4 read from file matrix4.txt
//chosen algorithm: LU pivoting 
int main(int argc, char* argv[])
{
    clock_t timer; 
    
    struct squareMatrix matrix4, inverse;
    FILE *f; 
        
    f = fopen("../Input/matrix4.txt", "r");
    readSquareMatrixFile(f, &matrix4);
    fclose(f);

    timer = clock();
    matrixInversePivoting(&matrix4, &inverse);
    timer = clock() - timer; 

    f = fopen("../Output/inverse.txt", "w");
    printSquareMatrixFile(f, &inverse);
    fclose(f);
    printf("Time to compute the inverse: %0.6f seconds\n", ((double)timer)/CLOCKS_PER_SEC);

    free(matrix4.mat);
    free(inverse.mat);
}
