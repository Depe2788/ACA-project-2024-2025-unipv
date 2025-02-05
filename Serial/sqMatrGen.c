#include "functions.h"

int main(){
    struct squareMatrix matrix4;
    FILE *f = fopen("../Input/matrix4.txt", "w");
    srand(time(NULL));
    matrix4.n = 1000; 
    matrix4.mat = (double *)malloc(matrix4.n * matrix4.n * sizeof(double));
    initializeMatrix(matrix4.mat, matrix4.n, matrix4.n);
    printSquareMatrixtFile(f, &matrix4);
    fclose(f);
    free(matrix4.mat);
}