#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int threshold = 1000;

void swap(int *a, int *b) {
    int tmp = *a;
    *a = *b;
    *b = tmp;
}

void partition(int *v, int* i, int* j, int low, int high) {
    *i = low;
    *j = high;
    int pivot = v[(low + high) / 2];
    do {
        while (v[*i] < pivot) (*i)++;
        while (v[*j] > pivot) (*j)--;
        if (i <= j) {
            swap(v + *i, v + *j);
            (*i)++;
            (*j)--;
        }
    } while (*i <= *j);
}

void quicksort_tasks(int *v, int low, int high) {
    int i, j;
    partition(v, &i, &j, low, high);
    if (high - low < threshold || (j - low < threshold || high - i < threshold)) {
        if (low < j)
            quicksort_tasks(v, low, j);
        if(i < high)
            quicksort_tasks(v, i, high);
    } else {
        #pragma omp task
            { quicksort_tasks(v, low, j); }
        quicksort_tasks(v, i, high);
    }
}

void initial (int *v, int size) {
    for (int i = 0; i < size; i++) {
        v[i] = rand() / RAND_MAX * 1000;
    }
}

int main() {
    for (threshold = 1000; threshold <= 100000; threshold *= 10) {
        char filename[100];
        sprintf(filename, "quik_sort_th=%d.csv", threshold);
        FILE *fout = fopen(filename, "w");
        for (int size = 1e6; size <= 1e8; size *= 10) {
            int *array = malloc(sizeof(int) * size);
            initial(array, size);
            double serial_time = -omp_get_wtime();
            quicksort_tasks(array, 0, size - 2);
            serial_time += omp_get_wtime();
            free(array);

            printf("n = %d\nSerial: %.6lf\n", size, serial_time);            

            for (int threads = 2; threads <= 8; threads += 2) {
                array = malloc(sizeof(int) * size);
                initial(array, size);

                double parallel_time = -omp_get_wtime();
                #pragma omp parallel num_threads(threads) 
                {
                    #pragma omp single 
                    quicksort_tasks(array, 0, size - 1);
                }
                parallel_time += omp_get_wtime();
                printf("N = %d, Threshold = %d, Threads = %d\nSpeedup: %.6lf\n", size, threshold, threads, serial_time / parallel_time);
                fprintf(fout, "%d; %d; %.4lf\n", size, threads, serial_time / parallel_time);
                free(array);
            }
        }
        fclose(fout);
    }
    return 0;
}
