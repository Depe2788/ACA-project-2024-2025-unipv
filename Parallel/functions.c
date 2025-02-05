#include "functions.h"

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

void readSquareMatrixFile(FILE *f, struct squareMatrix *matrix){
        if (f == NULL) {
                perror("Invalide file pointer");
                exit(1);
        }
        char buf[100]; 
        fgets(buf, sizeof(buf), f);
        sscanf(buf, "%i", &matrix->n);
        matrix->mat = (double *)malloc((matrix->n) * (matrix->n) * sizeof(double));
        for(int i = 0; i < ((matrix->n) * (matrix->n) ); i++){     
                fgets(buf, sizeof(buf),f);
                matrix->mat[i]=atof(buf);
        }
}

//one element on each row
void printSquareMatrixFile(FILE *f, struct squareMatrix *matrix)
{
        if (f == NULL) {
                perror("Invalide file pointer");
                exit(1);
        }
        fprintf(f, "%i\n", matrix->n);
        for(int i = 0; i < matrix->n; i++){   
                for(int j = 0; j < matrix->n; j++){   
                        fprintf(f, "%.6f\n", matrix->mat[matrix->n * i + j]);
                }
        }
}

void transposeInPlace(struct squareMatrix *matrix) {
    double temp;
    for (int i = 0; i < matrix->n; i++) {
        for (int j = i + 1; j < matrix->n; j++) {
            temp = matrix->mat[matrix->n * i + j];
            matrix->mat[matrix->n * i + j] = matrix->mat[matrix->n * j + i];
            matrix->mat[matrix->n * j + i] = temp;
        }
    }
}

void printIntegerVector(int *v, int n){
        for(int i = 0; i < n; i++){
                printf("%i ", v[i]);
        }
}

void printVectorStruct(struct vector *v){
        for(int i = 0; i < v->length; i++){
                printf("%f ", v->vect[i]);
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

//linear systems Ax = b; A must be a square matrix non singular
//forward substitution method for lower triangular systems 
//solution x will be in b
void forwardSubstitution(struct squareMatrix *A, struct vector *b){
        for (int i = 0; i < A->n; i++){
                for (int j = 0; j < i; j++){
                     b->vect[i] -= A->mat[(A->n) * i + j] * b->vect[j];
                }
                if (A->mat[(A->n) * i + i] == 0) {
                        printf("Error: diagonal element 0 in the forward substitution (A singular).\n");
                        exit(1);
                }
                b->vect[i] /= A->mat[(A->n) * i + i];
        }
}

//linear systems Ax = b; A must be a square matrix non singular
//backward substitution method for upper triangular systems 
void backwardSubstitution(struct squareMatrix *A, struct vector *b){
        for (int i = (A->n) - 1; i >= 0; i--){
                for (int j = i+1; j < A->n; j++){
                    b->vect[i] -= A->mat[(A->n) * i + j] * b->vect[j];
                }
                if (A->mat[(A->n) * i + i] == 0) {
                        printf("Error: diagonal element 0 in the backward substitution (A singular).\n");
                        exit(1);
                }
                b->vect[i] /= A->mat[(A->n) * i + i];
        }
}

//compute the inverse with LU pivoting
void matrixInversePivoting(struct squareMatrix *A, struct squareMatrix *inverse, int my_rank, int size, double * timer) {
        
        struct matrix inversePart;
        inversePart.ncols = A->n;

        int offset;

        double tmp;
        int maxIndex; 

        struct squareMatrix L, U, P; 
        L.n = A->n;
        U.n = A->n; 
        P.n = A->n; 
        L.mat = (double *)malloc(L.n * L.n * sizeof(double));
        U.mat = (double *)malloc(U.n * U.n * sizeof(double));
        P.mat = (double *)malloc(P.n * P.n * sizeof(double));

        struct vector pe;
        pe.length = A->n;
        pe.vect = (double *)malloc(pe.length * sizeof(double));

        int *recvcounts;
        int *displs;

        if (my_rank == 0) {
                //recvcounts[i] is the number of elements of inversePart of process i
                recvcounts = (int *)malloc(size * sizeof(int));
                displs = (int *)malloc(size * sizeof(int));

                int remaining = A->n % size;
                int div = A->n / size;
                for (int i = 0; i < size; i++) {
                        recvcounts[i] = div * A->n;
                        if (i < remaining) {
                                recvcounts[i] += A->n;  
                        }
                }

                displs[0] = 0;
                for (int i = 1; i < size; i++) {
                        displs[i] = displs[i - 1] + recvcounts[i - 1];
                }
                
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

                for (int k = 0; k < (A->n) - 1; k++){
                        tmp = fabs(U.mat[(U.n) * k + k]);
                        maxIndex = k;
                        for(int j = k + 1; j < (A->n); j++){
                                if(fabs(U.mat[(U.n) * j + k]) > tmp){
                                        tmp = fabs(U.mat[(U.n)*j+k]);
                                        maxIndex = j; 
                                }
                        }

                        // Controllo pivot nullo
                        if (tmp < 1e-10) {
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
                        /*
                        printf("matrix L iter k: %i:\n", k);
                        printMatrix(L.mat, L.n, L.n);
                        printf("matrix U iter k: %i:\n", k);
                        printMatrix(U.mat, U.n, U.n);
                        */
                }
        }

        //forward L, U and P to other processes
        MPI_Bcast(L.mat, L.n * L.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(U.mat, U.n * U.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(P.mat, P.n * P.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        /*
        if(my_rank == 0){
                printf("matrix L FINAL:\n");
                printMatrix(L.mat, L.n, L.n);
                printf("matrix U FINAL:\n");
                printMatrix(U.mat, U.n, U.n);
                printf("matrix P FINAL:\n");
                printMatrix(P.mat, P.n, P.n);
        }
        */

        MPI_Scatter(recvcounts, 1, MPI_INT, &inversePart.nrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
        inversePart.nrows /= inversePart.ncols; 
        inversePart.mat = (double *)malloc(inversePart.nrows * inversePart.ncols * sizeof(double));
        
        MPI_Scatter(displs, 1, MPI_INT, &offset, 1, MPI_INT, 0, MPI_COMM_WORLD);
        offset /= inversePart.ncols;
        //offset is the index of the initial row of InversePart with respect to the inverse
        //printf("process %d, offset = %d\n", my_rank, offset);
       
        for (int i = 0; i < inversePart.nrows; i++){
                for(int k = 0; k < pe.length; k++){
                        pe.vect[k] = P.mat[P.n * k + i + offset];
                }
                forwardSubstitution(&L, &pe);
                backwardSubstitution(&U, &pe);
                //copy pe (ci) in the inversePart
                for(int k = 0; k < inversePart.ncols; k++){
                        inversePart.mat[inversePart.ncols * i + k] = pe.vect[k];         
                }
        }

        MPI_Gatherv(inversePart.mat, inversePart.nrows * inversePart.ncols, MPI_DOUBLE, inverse->mat, 
            recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if(my_rank == 0){
                transposeInPlace(inverse);
                free(recvcounts);
                free(displs);
        }

        free(L.mat);
        free(U.mat);
        free(P.mat);
        free(pe.vect);
        free(inversePart.mat);
}


//compute the inverse with LU pivoting
void matrixInversePivotingImproved(struct squareMatrix *A, struct squareMatrix *inverse, int my_rank, int size, double * timer) {
        
        struct matrix inversePart;
        inversePart.ncols = A->n;
        
        int offset;

        double tmp;
        int maxIndex; 

        int *recvcounts;
        int *displs;
        int *sendcountsl;
        int *displsl;

        struct matrix lines;
        lines.ncols = A->n;

        struct vector lPart;

        int rowsbasic, remaining, div; 

        struct squareMatrix L, U, P; 
        L.n = A->n;
        U.n = A->n; 
        P.n = A->n; 
        L.mat = (double *)malloc(L.n * L.n * sizeof(double));
        U.mat = (double *)malloc(U.n * U.n * sizeof(double));
        P.mat = (double *)malloc(P.n * P.n * sizeof(double));

        struct vector pe;
        pe.length = A->n;
        pe.vect = (double *)malloc(pe.length * sizeof(double));
        
        if(my_rank == 0){

                recvcounts = (int *)malloc(size * sizeof(int));
                displs = (int *)malloc(size * sizeof(int));

                sendcountsl = (int *)malloc(size * sizeof(int));
                displsl = (int *)malloc(size * sizeof(int));

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
        }

        for (int k = 0; k < (A->n) - 1; k++){    

                if (my_rank == 0){
                        tmp = fabs(U.mat[(U.n) * k + k]);
                        maxIndex = k;
                        for(int j = k + 1; j < (A->n); j++){
                                if(fabs(U.mat[(U.n) * j + k]) > tmp){
                                        tmp = fabs(U.mat[(U.n)*j+k]);
                                        maxIndex = j; 
                                }
                        }

                        // Controllo pivot nullo
                        if (tmp < 1e-10) {
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
                        //L is transpose so take it into account (s and maxIndex are inverted)
                        if(k >= 1){
                                for(int s = 0; s < k; s++){
                                        tmp = L.mat[(L.n) * s + k];
                                        L.mat[(L.n) * s + k] = L.mat[(L.n) * s + maxIndex];
                                        L.mat[(L.n) * s + maxIndex] = tmp;
                                }
                                
                        }
                        
                        rowsbasic = (U.n - k - 1) / size;
                        remaining = (U.n - k - 1) % size;
                        for (int i = 0; i < size; i++) {
                                sendcountsl[i] = rowsbasic;
                                if (i < remaining) {
                                        sendcountsl[i] += 1;  
                                }
                                recvcounts[i] = sendcountsl[i] * lines.ncols;
                        }

                        displs[0] = 0;
                        displsl[0] = 0;
                        for (int i = 1; i < size; i++) {
                                displs[i] = displs[i - 1] + recvcounts[i - 1];
                                displsl[i] = displsl[i - 1] + sendcountsl[i - 1];
                        }
                        /*
                        printf("!!!!!matrix U k=%i after changing rows:\n", k);
                        printMatrix(U.mat, U.n, U.n);
                        printf("\n\n\n");
                        */
                }

                MPI_Scatter(sendcountsl, 1, MPI_INT, &lines.nrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
                lines.mat = (double *)malloc(lines.nrows * lines.ncols * sizeof(double));
                lPart.length = lines.nrows;
                lPart.vect = (double *)malloc(lPart.length * sizeof(double));

                //Process 0 sends a part of U under the pivot line to each process (lines.mat)
                MPI_Scatterv(&U.mat[U.n * (k + 1)], recvcounts, displs, MPI_DOUBLE, lines.mat, lines.nrows * lines.ncols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                //Process 0 sends the pivot line in broadcast
                MPI_Bcast(&U.mat[U.n * k + k], lines.ncols, MPI_DOUBLE, 0, MPI_COMM_WORLD);


                for (int i = 0; i < lines.nrows; i++){
                        lPart.vect[i] = lines.mat[lines.ncols * i + k] / U.mat[U.n * k + k];
                        //lines.mat[lines.ncols * i] = 0;
                        for (int s = 0; s < lines.ncols - k; s++){
                                lines.mat[lines.ncols * i + k + s] = lines.mat[lines.ncols * i + k + s] - lPart.vect[i] * U.mat[U.n * k + k + s];
                        }
                }
                
                //printf("vector lPart!!\n");
                //printVectorStruct(&lPart);
                //printf("\n");

                //rimettere nella U 
                MPI_Gatherv(lines.mat, lines.nrows * lines.ncols, MPI_DOUBLE, &U.mat[U.n * (k + 1)], recvcounts, displs,  MPI_DOUBLE, 0, MPI_COMM_WORLD);
                //metto lPart nella L 
                MPI_Gatherv(lPart.vect, lPart.length, MPI_DOUBLE, &L.mat[L.n * k + (k+1)], sendcountsl, displsl,  MPI_DOUBLE, 0, MPI_COMM_WORLD);
                /*
                if(my_rank == 0){
                        printf("matrix L k=%i:\n", k);
                        printMatrix(L.mat, L.n, L.n);
                        printf("matrix U k=%i:\n", k);
                        printMatrix(U.mat, U.n, U.n);
                        
                }
                */
                free(lines.mat);
                free(lPart.vect);
        }

        if(my_rank == 0){
                transposeInPlace(&L);

                free(sendcountsl);
                free(displsl);
        }

        //forward L, U and P to other processes
        MPI_Bcast(L.mat, L.n * L.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(U.mat, U.n * U.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(P.mat, P.n * P.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        /*
        if(my_rank == 0){
                printf("matrix L FINAL:\n");
                printMatrix(L.mat, L.n, L.n);
                printf("matrix U FINAL:\n");
                printMatrix(U.mat, U.n, U.n);
                printf("matrix P FINAL:\n");
                printMatrix(P.mat, P.n, P.n);
        }
        */

       //--------
        if (my_rank == 0) {

                remaining = A->n % size;
                div = A->n / size;
                for (int i = 0; i < size; i++) {
                        recvcounts[i] = div * A->n;
                        if (i < remaining) {
                                recvcounts[i] += A->n;  
                        }
                }

                displs[0] = 0;
                for (int i = 1; i < size; i++) {
                        displs[i] = displs[i - 1] + recvcounts[i - 1];
                }
        }

        MPI_Scatter(recvcounts, 1, MPI_INT, &inversePart.nrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
        inversePart.nrows /= inversePart.ncols; 
        inversePart.mat = (double *)malloc(inversePart.nrows * inversePart.ncols * sizeof(double));
        
        MPI_Scatter(displs, 1, MPI_INT, &offset, 1, MPI_INT, 0, MPI_COMM_WORLD);
        offset /= inversePart.ncols;

        //------------------------------

        for (int i = 0; i < inversePart.nrows; i++){
                for(int k = 0; k < pe.length; k++){
                        pe.vect[k] = P.mat[P.n * k + i + offset];
                }
                forwardSubstitution(&L, &pe);
                backwardSubstitution(&U, &pe);
                //copy pe (ci) in the inversePart
                for(int k = 0; k < inversePart.ncols; k++){
                        inversePart.mat[inversePart.ncols * i + k] = pe.vect[k];         
                }
        }

        MPI_Gatherv(inversePart.mat, inversePart.nrows * inversePart.ncols, MPI_DOUBLE, inverse->mat, 
            recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if(my_rank == 0){
                transposeInPlace(inverse);
                free(recvcounts);
                free(displs);
        }

        free(L.mat);
        free(U.mat);
        free(P.mat);
        free(pe.vect);
        free(inversePart.mat);
}
