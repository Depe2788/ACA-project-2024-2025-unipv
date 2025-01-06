#include "functions.h"

//Inside the function there is the dynamic allocation, so the address mat changes, pass the argument by reference
void readMatrixFile(FILE *f, double **mat, int *nrows, int *ncols){
        if (f == NULL) {
                perror("Invalide file pointer");
                exit(1);
        }
        char buf[100]; 
        fgets(buf, sizeof(buf), f);
        sscanf(buf, "%i %i", nrows, ncols);
        *mat = (double *)malloc((*nrows) * (*ncols) * sizeof(double));
        for(int i = 0; i < ((*nrows) * (*ncols)); i++){     
                fgets(buf, sizeof(buf),f);
                (*mat)[i]=atof(buf);
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

void printMatrix(double *mat, int nrows, int ncols)
{
        for (int i = 0; i < nrows; i++){
                for(int j = 0; j < ncols; j++){
                        printf("%0.2f ", mat[ncols*i+j]);
                } 
                printf("\n");
        }
}

void printMatrixFile(FILE *f, double *mat, int nrows, int ncols)
{
        if (f == NULL) {
                perror("Invalide file pointer");
                exit(1);
        }
        fprintf(f, "%i %i\n", nrows, ncols);
        for(int i = 0; i < nrows; i++){   
                for(int j = 0; j < ncols; j++){   
                        fprintf(f, "%.2f ", mat[ncols * i + j]);
                }
                fprintf(f, "\n");
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

//linear systems Ax = b; A must be a square matrix non singular
//forward substitution method for lower triangular systems 
void forwardSubstitution(struct squareMatrix *A, struct vector *b, struct vector *x){
        for (int i = 0; i < A->n; i++){
                x->vect[i] = b->vect[i];
                for (int j = 0; j < i; j++){
                     x->vect[i] -= A->mat[(A->n) * i + j] * x->vect[j];
                }
                if (A->mat[(A->n) * i + i] == 0) {
                        printf("Error: diagonal element 0 in the forward substitution (A singular).\n");
                        exit(1);
                }
                x->vect[i] /= A->mat[(A->n) * i + i];
        }
}

//linear systems Ax = b; A must be a square matrix non singular
//backward substitution method for upper triangular systems 
void backwardSubstitution(struct squareMatrix *A, struct vector *b, struct vector *x){
        for (int i = (A->n) - 1; i >= 0; i--){
                x->vect[i] = b->vect[i];
                for (int j = i+1; j < A->n; j++){
                    x->vect[i] -= A->mat[(A->n) * i + j] * x->vect[j];
                }
                if (A->mat[(A->n) * i + i] == 0) {
                        printf("Error: diagonal element 0 in the backward substitution (A singular).\n");
                        exit(1);
                }
                x->vect[i] /= A->mat[(A->n) * i + i];
        }
}

//compute the inverse with LU pivoting
void matrixInversePivoting(struct squareMatrix *A, struct squareMatrix *inverse){
        inverse->n = A->n;
        inverse->mat = (double *)malloc(inverse->n * inverse->n * sizeof(double));

        //PA=LU find L, U, P
        struct squareMatrix L, U, P;
        L.n = A->n;
        L.mat = (double *)malloc(L.n * L.n * sizeof(double));
        U.n = A->n; 
        U.mat = (double *)malloc(U.n * U.n * sizeof(double));
        P.n = A->n; 
        P.mat = (double *)malloc(P.n * P.n * sizeof(double));

        //initialize L and P (identity matrices) and U = A
        for (int i = 0; i < A->n; i++){
                for(int j = 0; j < A->n; j++){  
                        if(i == j){
                                L.mat[(L.n) * i + j] = 1;
                                P.mat[(P.n) * i + j] = 1;
                        } else {
                                L.mat[(L.n) * i + j] = 0;
                                P.mat[(P.n) * i + j] = 0;
                        } 
                        U.mat[(U.n) * i + j] = A->mat[(A->n) * i + j];
                }
        }

        double tmp;
        int maxIndex; 
        for (int k = 0; k < (A->n) - 1; k++){
                tmp = U.mat[(U.n)*k+k];
                maxIndex = k;
                for(int j = k; j < (A->n); j++){
                        if(fabs(U.mat[(U.n) * j + k]) > fabs(tmp)){
                                tmp = U.mat[(U.n)*j+k];
                                maxIndex = j; 
                        }
                }

                //check if the pivot is equal to 0
                if (fabs(tmp) < 1e-10) {
                        printf("Error: Pivot 0 or too small.\n");
                        exit(EXIT_FAILURE);
                }

                for(int s = k; s < A->n; s++){
                        tmp = U.mat[(U.n) * k + s];
                        U.mat[(U.n) * k + s] = U.mat[(U.n) * maxIndex + s]; 
                        U.mat[(U.n) * maxIndex + s] = tmp;
                }

                for(int s = 0; s < A->n; s++){
                        tmp = P.mat[(P.n) * k + s];
                        P.mat[(P.n) * k + s] = P.mat[(P.n) * maxIndex + s]; 
                        P.mat[(P.n) * maxIndex + s] = tmp;
                }

                if(k >= 1){
                        for(int s = 0; s < k; s++){
                                tmp = L.mat[(L.n) * k + s];
                                L.mat[(L.n) * k + s] = L.mat[(L.n) * maxIndex + s];
                                L.mat[(L.n) * maxIndex + s] = tmp;
                        }
                        
                }

                for (int i = k+1; i < A->n; i++){
                        L.mat[L.n * i + k] = U.mat[U.n * i + k] / U.mat[U.n * k + k];
                        for (int s = k; s < A->n; s++){
                                U.mat[U.n * i + s] = U.mat[U.n * i + s] - L.mat[L.n * i + k] * U.mat[U.n * k + s];
                        }
                }
        }  
        
        //output: L, U, P
        /*
        printf("Matrix L:\n");
        printMatrix(L.mat, L.n, L.n);
        printf("Matrix U:\n");
        printMatrix(U.mat, U.n, U.n); 
        printf("Matrix P:\n");
        printMatrix(P.mat, P.n, P.n); 
        */
        
        struct vector y; 
        y.length = A->n;
        y.vect = (double *)malloc(y.length * sizeof(double));
        struct vector pe;
        pe.length = A->n;
        pe.vect = (double *)malloc(pe.length * sizeof(double));
        struct vector c;
        c.length = A->n;
        c.vect = (double *)malloc(c.length * sizeof(double));

        for (int i = 0; i < A->n; i++){
                for(int k = 0; k < pe.length; k++){
                        pe.vect[k] = P.mat[P.n * k + i];
                }
                forwardSubstitution(&L, &pe, &y);
                backwardSubstitution(&U, &y, &c);
                //copy c in the inverse
                for(int k = 0; k < A->n; k++){
                        inverse->mat[(A->n) * k + i] = c.vect[k];
                }
        }
  
        free(L.mat);
        free(U.mat);
        free(P.mat);
        free(y.vect);
        free(pe.vect);
        free(c.vect);
}
