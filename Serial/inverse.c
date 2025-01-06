#include "functions.h"

//compute the inverse of a square matrix matrix4 read from file matrix4.txt
int main(int argc, char* argv[])
{
    clock_t timer; 
    
    struct squareMatrix matrix4, inverse;
    FILE *f; 
        
    f = fopen("matrix4.txt", "r");
    readMatrixFile(f, &matrix4.mat, &matrix4.n, &matrix4.n);
    fclose(f);

    timer = clock();
    matrixInversePivoting(&matrix4, &inverse);
    timer = clock() - timer; 

    f = fopen("inverse.txt", "w");
    printMatrixFile(f, inverse.mat, inverse.n, inverse.n);
    fclose(f);
    printf("Time to compute the inverse: %0.6f seconds\n", ((double)timer)/CLOCKS_PER_SEC);
}
