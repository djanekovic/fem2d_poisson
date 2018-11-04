#include <stdbool.h>
#include <stdio.h>

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

void print_global_RHS(REAL *matrix, int num_elements)
{
    for (int i = 0; i < num_elements; i++) {
        printf ("%lf\n", matrix[i]);
    }
}

void print_Ax_b(REAL *A, REAL *b, int num_elements)
{
    for (int i = 0; i < num_elements; i++) {
        for (int j = 0; j < num_elements; j++) {
            printf ("%lf  ", A[i*num_elements + j]);
        }
        if (i == num_elements/2)
            printf ("=");
        else
            printf (" ");
        printf("  %lf\n", b[i]);
    }
}
/*
void report_triangulateio(struct triangulateio *out,
                          bool points,
                          bool triangles)
{
    if (points) {
        for (size_t i = 0; i < out->numberofpoints; i++) {
            printf("Point %d [%lf, %lf] - %d\n",
                   i,
                   out->pointlist[2 * i],
                   out->pointlist[2 * i + 1],
                   out->pointmarkerlist[i]);
        }
    }

    if (triangles) {
        for (size_t i = 0; i < out->numberoftriangles; i++) {
            unsigned int point_id[3];
            memcpy(point_id, out->trianglelist[i * 3], 3 * sizeof(int));
            printf("Triangle %d:\n", i);
            printf("\tPoints %d [%lf, %lf], %d [%lf %lf], %d [%lf %lf]\n",
                   point_id[0],
                   out->pointlist[point_id[0]],
                   out->pointlist[point_id[0] + 1],
                   point_id[1],
                   out->pointlist[point_id[1]],
                   out->pointlist[point_id[1] + 1],
                   point_id[2],
                   out->pointlist[point_id[2]],
                   out->pointlist[point_id[2] + 1]);
        }
    }
}
*/
