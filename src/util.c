#include <stdio.h>
#include <stdbool.h>

#include "util.h"

void print_local_matrix(REAL *matrix)
{
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%lf ", matrix[i * 3 + j]);
        }
        printf("\n");
    }
}

void print_global_matrix(REAL *matrix, int num_elements, bool num)
{
    for (int i = 0; i < num_elements; i++) {
        for (int j = 0; j < num_elements; j++) {
            if (num)
                printf("%lf ", matrix[i * num_elements + j]);
            else
                printf("%c ", (matrix[i * num_elements + j]) ? 'X' : '0');
        }
        printf("\n");
    }
}
