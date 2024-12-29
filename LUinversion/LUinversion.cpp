#include <stdio.h>
#include <stdlib.h>
#include <time.h>  // Libreria per misurare il tempo di esecuzione

// Funzione per leggere una matrice quadrata da un file
void readMatrix(double** matrix, FILE* file, int n) {

    int i, j;

    for (i = 0; i < n; i++) {
        //Alloco spazio per i singoli elementi
        matrix[i] = (double*)malloc(n * sizeof(double));
        if (matrix[i] == NULL) {
            printf("Errore: impossibile allocare memoria per la riga %d.\n", i);
            exit(EXIT_FAILURE);
        }
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (fscanf_s(file, "%lf", &matrix[i][j]) != 1) {
                printf("Errore: lettura del file fallita alla posizione (%d, %d).\n", i, j);
                exit(EXIT_FAILURE);
            }
        }
    }
}

// Funzione per memorizzare una matrice in un file, elemento per elemento 
void printMatrix(double** matrix, int n, FILE* file) {
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            fprintf(file, "%lf ", matrix[i][j]);
        }
        fprintf(file, "\n");
    }
}

// Funzione per eseguire la fattorizzazione LU
void luDecomposition(double** matrix, double** L, double** U, int size) {
    clock_t start = clock(); 

    //Outer loop
    for (int i = 0; i < size; i++) {

        //First inner loop
        for (int j = i; j < size; j++) {
            U[i][j] = matrix[i][j];         
            for (int k = 0; k < i; k++) {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }

        //Second inner loop
        for (int j = i; j < size; j++) {
            if (i == j)
                L[i][i] = 1;
            else {
                L[j][i] = matrix[j][i];
                for (int k = 0; k < i; k++) {
                    L[j][i] -= L[j][k] * U[k][i];
                }
                L[j][i] /= U[i][i];
            }
        }
    }

    clock_t end = clock(); 
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Decomposizione LU: %f secondi\n", time_taken);
}

// Funzione per risolvere il sistema LY = B (sostituzione in avanti)
void forwardSubstitution(double** L, double* B, double* Y, int size) {
    for (int i = 0; i < size; i++) {
        Y[i] = B[i];
        for (int j = 0; j < i; j++) {
            Y[i] -= L[i][j] * Y[j];
        }
    }
}

// Funzione per risolvere il sistema UX = Y (sostituzione all'indietro)
void backwardSubstitution(double** U, double* Y, double* X, int size) {
    for (int i = size - 1; i >= 0; i--) {
        X[i] = Y[i];
        for (int j = i + 1; j < size; j++) {
            X[i] -= U[i][j] * X[j];
        }
        X[i] /= U[i][i];
    }
}

// Funzione per calcolare l'inversa della matrice
void invertMatrixLU(double** matrix, double** inverse, int size) {
    double** L = (double**)malloc(size * sizeof(double*));
    double** U = (double**)malloc(size * sizeof(double*));
    for (int i = 0; i < size; i++) {
        L[i] = (double*)calloc(size, sizeof(double));
        U[i] = (double*)calloc(size, sizeof(double));
    }

    luDecomposition(matrix, L, U, size);

    double* B = (double*)malloc(size * sizeof(double));
    double* Y = (double*)malloc(size * sizeof(double));
    double* X = (double*)malloc(size * sizeof(double));

    clock_t start = clock();


    for (int col = 0; col < size; col++) {
        for (int i = 0; i < size; i++) {
            B[i] = (i == col) ? 1.0 : 0.0;
        }

        forwardSubstitution(L, B, Y, size);
        backwardSubstitution(U, Y, X, size);

        for (int i = 0; i < size; i++) {
            inverse[i][col] = X[i];
        }
    }

    clock_t end = clock();  // Fine del profiling
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Forward e backward: %f secondi\n", time_taken);


    for (int i = 0; i < size; i++) {
        free(L[i]);
        free(U[i]);
    }
    free(L);
    free(U);
    free(B);
    free(Y);
    free(X);
}

int main(int argc, char* argv[]) {

    if (argc != 2) {
        exit(EXIT_FAILURE);
    }

    FILE* file;
    if (fopen_s(&file, argv[1], "r") != 0) {
        printf("Errore: impossibile aprire il file '%s'.\n", argv[1]);
        exit(EXIT_FAILURE);
    }

    int n;

    if (fscanf_s(file, "%d", &n) != 1) {
        printf("Errore: formato del file non valido (manca la dimensione della matrice).\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    // Allocazione della matrice. 
    double** matrix = (double**)malloc(n * sizeof(double*));
    if (matrix == NULL) {
        printf("Errore: memoria insufficiente per la matrice.\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    // Lettura della matrice dal file
    readMatrix(matrix, file, n);
    fclose(file);

    // Allocazione della matrice inversa
    // Alloco spazio per un vettore di n puntatori.
    double** inverse = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        //Ogni cella del vettore contiene un puntatore che punta ad una riga della matrice
        inverse[i] = (double*)malloc(n * sizeof(double));
    }

    // Inizio della misurazione del tempo
    clock_t start = clock();

    // Calcolo dell'inversa
    invertMatrixLU(matrix, inverse, n);

    // Fine della misurazione del tempo
    clock_t end = clock();

    // Calcolo il tempo di esecuzione in secondi
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;

    printf("Tempo di esecuzione per l'inversione della matrice: %f secondi\n", time_taken);

    // Scrittura della matrice inversa nel file
    FILE* resultFile;
    if (fopen_s(&resultFile, "result.txt", "w") != 0) {
        perror("Errore nell'aprire il file result.txt");
        exit(EXIT_FAILURE);
    }

    printMatrix(inverse, n, resultFile);
    fclose(resultFile);

    // Deallocazione della memoria
    for (int i = 0; i < n; i++) {
        free(matrix[i]);
        free(inverse[i]);
    }
    free(matrix);
    free(inverse);

    return 0;
}
