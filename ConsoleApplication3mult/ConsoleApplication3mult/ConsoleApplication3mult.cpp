#include <stdio.h>
#include <stdlib.h>
#include <time.h>  // Aggiungi questa libreria per misurare il tempo

// Funzione per leggere una matrice quadrata da un file
void readMatrix(double** matrix, FILE* file, int n) {
    int i, j;

    for (i = 0; i < n; i++) {
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

// Funzione per stampare una matrice su un file
void printMatrix(double** matrix, int n, FILE* file) {
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            fprintf(file, "%lf ", matrix[i][j]);
        }
        fprintf(file, "\n");
    }
}

// Funzione per eseguire la moltiplicazione tra due matrici quadrate
void matrixMul(double** m1, double** m2, double** m3, int n) {
    int i, j, k;

    // Inizializzo la matrice risultante a zero
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            m3[i][j] = 0.0;
        }
    }

    // Moltiplicazione delle matrici
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                m3[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        printf("Uso: %s <file_matrice1> <file_matrice2>\n", argv[1]);
        printf("Uso: %s <file_matrice1> <file_matrice2>\n", argv[2]);
        exit(EXIT_FAILURE);
    }

    FILE* file1, * file2;
    if (fopen_s(&file1, argv[1], "r") != 0) {
        printf("Errore: impossibile aprire il file '%s'.\n", argv[1]);
        exit(EXIT_FAILURE);
    }

    if (fopen_s(&file2, argv[2], "r") != 0) {
        printf("Errore: impossibile aprire il file '%s'.\n", argv[2]);
        exit(EXIT_FAILURE);
    }

    int n;

    // Leggi la dimensione delle matrici
    if (fscanf_s(file1, "%d", &n) != 1 || fscanf_s(file2, "%d", &n) != 1) {
        printf("Errore: formato del file non valido (manca la dimensione della matrice).\n");
        fclose(file1);
        fclose(file2);
        exit(EXIT_FAILURE);
    }

    // Allocazione delle matrici
    double** matrix1 = (double**)malloc(n * sizeof(double*));
    double** matrix2 = (double**)malloc(n * sizeof(double*));
    double** result = (double**)malloc(n * sizeof(double*));

    for (int i = 0; i < n; i++) {
        matrix1[i] = (double*)malloc(n * sizeof(double));
        matrix2[i] = (double*)malloc(n * sizeof(double));
        result[i] = (double*)malloc(n * sizeof(double));
    }

    // Lettura delle matrici dai file
    readMatrix(matrix1, file1, n);
    readMatrix(matrix2, file2, n);

    fclose(file1);
    fclose(file2);

    // Misurazione del tempo di esecuzione
    clock_t start = clock();

    // Esegui la moltiplicazione delle matrici
    matrixMul(matrix1, matrix2, result, n);

    clock_t end = clock();

    // Calcolo il tempo di esecuzione in secondi
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tempo di esecuzione per la moltiplicazione della matrice: %f secondi\n", time_taken);

    // Scrittura del risultato in un file
    FILE* resultFile;
    if (fopen_s(&resultFile, "result.txt", "w") != 0) {
        printf("Errore nell'aprire il file result.txt.\n");
        exit(EXIT_FAILURE);
    }

    printMatrix(result, n, resultFile);
    fclose(resultFile);

    // Deallocazione della memoria
    for (int i = 0; i < n; i++) {
        free(matrix1[i]);
        free(matrix2[i]);
        free(result[i]);
    }
    free(matrix1);
    free(matrix2);
    free(result);

    return 0;
}
