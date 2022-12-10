//mpiexec -n 4 lab1.exe 
#include <stdio.h>
#include <math.h>
#include<time.h>
#include<stdlib.h>

#include "mpi.h"

#define N 240

void print_results(const char prompt[5], unsigned long int a[N][N])
{
    int i, j;

    printf("\n%s\n", prompt);
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf(" %d", a[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

int main(int argc, char** argv)
{

    /* Инициализация библиотеки MPI */
    MPI_Init(&argc, &argv);

    double time_start, time_finish;
    time_start = MPI_Wtime();

    int i, j, rank, size, sum = 0;
    unsigned long int c[N][N];
    unsigned long int aa[N], cc[N];

    unsigned long int a[N][N];
    unsigned long int b[N][N];

    for ( i = 0; i < N; i++)
    {
        for ( j = 0; j < N; j++)
        {
            a[i][j] = rand() % 11;
            b[i][j] = rand() % 11;
        }
    }


    // Получить размер коммуникатора MPI_COMM_WORLD
    // (общее число процессов в рамках задачи)
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // Получить номер текущего процесса в рамках 
    // коммуникатора MPI_COMM_WORLD
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //транслировать вторую матрицу всем процессам
    MPI_Bcast(b, N * N, MPI_LONG, 0, MPI_COMM_WORLD);
    //распределить строки первой матрицы по разным процессам     
    if (size == N) {
        MPI_Scatter(a, N, MPI_UNSIGNED_LONG, aa, N , MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        

        //синхронизация процессов
        MPI_Barrier(MPI_COMM_WORLD);

        //умножение всеми процессами
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                sum = sum + aa[j] * b[j][i];
            }
            cc[i] = sum;
            sum = 0;
        }
        //собрать элементы в корневой процесс
        MPI_Gather(cc, N, MPI_UNSIGNED_LONG, c, N, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

    }
    else {
        int obhod = 0;
        unsigned long int prom_rez[N][N];
        if (size < N) {
            unsigned long int* cur_mass = (unsigned long int*)malloc(size * N * sizeof(unsigned long int));
            while (obhod < N / size) {
                for (i = 0; i < size; i++) {
                    for (j = 0; j < N; j++) {
                        cur_mass[i * N + j] = a[obhod * size + i][j];
                    }
                }
                MPI_Scatter(cur_mass, N, MPI_UNSIGNED_LONG, aa, N, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
                //синхронизация процессов
                MPI_Barrier(MPI_COMM_WORLD);
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        sum = sum + aa[j] * b[j][i];
                    }
                    cc[i] = sum;
                    sum = 0;
               
                }
                //собрать элементы в корневой процесс
                MPI_Gather(cc, N, MPI_UNSIGNED_LONG, prom_rez, N, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
                for (i = 0; i < size; i++) {
                    for (j = 0; j < N; j++) {
                        c[obhod * size + i][j] = prom_rez[i][j];
                    }
                }
                MPI_Barrier(MPI_COMM_WORLD);
                obhod += 1;
            }
            free(cur_mass);
        }
        else {
 
        }
    }

 
    
    if (rank == 0) {
        //раскомментировать, чтобы вывести матрицы
        /*print_results("A: ", a);
        print_results("B: ", b);
        print_results("C: ", c);*/
        time_finish = MPI_Wtime();
        printf("process: %d, time: %f\n", rank, time_finish - time_start);
    }
    
    // Освобождение подсистемы MPI
    MPI_Finalize();
    return 0;
}