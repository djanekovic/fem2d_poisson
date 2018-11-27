#ifndef UTIL_H

#define REAL double

#define DEBUG(global_A, global_F, num_of_points)                    \
    do {                                                            \
        printf ("\n");                                              \
        print_Ax_b(global_A, global_F, num_of_points);              \
    } while (0)


void print_local_matrix(REAL *matrix);
void print_global_matrix(REAL *matrix, int num_elements, bool num);
void print_global_RHS(REAL *matrix, int num_elements);
void print_Ax_b(REAL *A, REAL *b, int num_elements);

#endif /* UTIL_H */
