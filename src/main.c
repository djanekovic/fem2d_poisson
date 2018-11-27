/**
 * Finite element method solver. This solver solves laplace(u) = 6 on square
 * [0, 1]x[0, 1]. Boundary conditions are set to be 1 + x^2 + 2y^2.
 *
 * Exact solution of this problem is 1 + x^2 + 2y^2. We will test against this
 * solution in the end.
 *
 * Weak formulation is: ...
 *
 * Should take argument:
 *  - about domain (maybe file)
 *  - arguments that Triangle implements.
 *  - sparse implementation
 *  - solver choice
 *  - debug
 */

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lapacke.h>

#include "util.h"
#include "triangle/triangle.h"

#define F -6.0

#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)


#define COMPUTE_LOCAL_A(point_id, pointlist, matrix)                \
    do {                                                            \
        REAL x1 = pointlist[2 * point_id[0]];                       \
        REAL y1 = pointlist[2 * point_id[0] + 1];                   \
        REAL x2 = pointlist[2 * point_id[1]];                       \
        REAL y2 = pointlist[2 * point_id[1] + 1];                   \
        REAL x3 = pointlist[2 * point_id[2]];                       \
        REAL y3 = pointlist[2 * point_id[2] + 1];                   \
        REAL dx23 = x2 - x3;                                        \
        REAL dy23 = y2 - y3;                                        \
        REAL dx31 = x3 - x1;                                        \
        REAL dy31 = y3 - y1;                                        \
        REAL dx12 = x1 - x2;                                        \
        REAL dy12 = y1 - y2;                                        \
        REAL area = 0.5 * (dx31 * dy12 - dy31 * dx12);              \
        matrix[0] = 0.25 * (dx23 * dx23 + dy23 * dy23) / area;      \
        matrix[1] = 0.25 * (dx23 * dx31 + dy23 * dy31) / area;      \
        matrix[2] = 0.25 * (dx23 * dx12 + dy23 * dy12) / area;      \
        matrix[3] = 0.25 * (dx31 * dx23 + dy31 * dy23) / area;      \
        matrix[4] = 0.25 * (dx31 * dx31 + dy31 * dy31) / area;      \
        matrix[5] = 0.25 * (dx31 * dx12 + dy31 * dy12) / area;      \
        matrix[6] = 0.25 * (dx12 * dx23 + dy12 * dy23) / area;      \
        matrix[7] = 0.25 * (dx12 * dx31 + dy12 * dy31) / area;      \
        matrix[8] = 0.25 * (dx12 * dx12 + dy12 * dy12) / area;      \
    } while(0)

#define COMPUTE_LOCAL_F(point_id, pointlist, matrix)                \
    do {                                                            \
        REAL x1 = pointlist[2 * point_id[0]];                       \
        REAL y1 = pointlist[2 * point_id[0] + 1];                   \
        REAL x2 = pointlist[2 * point_id[1]];                       \
        REAL y2 = pointlist[2 * point_id[1] + 1];                   \
        REAL x3 = pointlist[2 * point_id[2]];                       \
        REAL y3 = pointlist[2 * point_id[2] + 1];                   \
        REAL dx31 = x3 - x1;                                        \
        REAL dy31 = y3 - y1;                                        \
        REAL dx12 = x1 - x2;                                        \
        REAL dy12 = y1 - y2;                                        \
        REAL area = 0.5 * (dx31 * dy12 - dy31 * dx12);              \
        REAL c_diag = area / 6.0;                                   \
        REAL c_off = area / 12.0;                                   \
        matrix[0] = c_diag; matrix[1] = c_off;  matrix[2] = c_off;  \
        matrix[3] = c_off;  matrix[4] = c_diag; matrix[5] = c_off;  \
        matrix[6] = c_off;  matrix[7] = c_off;  matrix[8] = c_diag; \
    } while(0)

static void generate_mesh(struct triangulateio *io, struct triangulateio *out);
static REAL boundary_condition(uint point_id, struct triangulateio *out);
static REAL exact_solution(uint point_id, struct triangulateio *out);

int main(int argc, char **argv)
{
    struct triangulateio in, out;
    uint num_of_points;
    uint num_of_triangles;
    REAL *local_A, *local_F, *global_A, *global_F;

    /* generate geometry */
    generate_mesh(&in, &out);

    double tmp[] = {
        0.0, 0.0,
        0.5, 0.0,
        1.0, 0.0,
        0.0, 0.5,
        0.5, 0.5,
        1.0, 0.5,
        0.0, 1.0,
        0.5, 1.0,
        1.0, 1.0
    };
    //memcpy(out.pointlist, tmp, sizeof(tmp));

    int tmp2[] = {
        0, 1, 3,
        4, 3, 1,
        1, 2, 4,
        5, 4, 2,
        3, 4, 6,
        7, 6, 4,
        4, 5, 7,
        8, 7, 5};
    //memcpy(out.trianglelist, tmp2, sizeof(tmp2));

    num_of_points = out.numberofpoints;
    num_of_triangles = out.numberoftriangles;

    /**
     * Allocate all matrix data structures
     */
    global_A = malloc(num_of_points * num_of_points * sizeof(REAL));
    global_F = malloc(num_of_points * sizeof(REAL));
    local_A = malloc(3 * 3 * sizeof(REAL));
    local_F = malloc(3 * 3 * sizeof(REAL));

    assert(global_A != NULL);
    assert(global_F != NULL);
    assert(local_A != NULL);
    assert(local_F != NULL);

    memset(global_A, 0, num_of_points * num_of_points * sizeof(REAL));
    memset(global_F, 0, num_of_points * sizeof(REAL));

    /**
     * Matrix assembly without boundary conditions
     */
    for (size_t n = 0; n < num_of_triangles; n++) {
        // three points that describe triangle
        uint point_id[3];
        memcpy(point_id, &out.trianglelist[n * 3], 3 * sizeof(int));
        COMPUTE_LOCAL_A(point_id, out.pointlist, local_A);
        COMPUTE_LOCAL_F(point_id, out.pointlist, local_F);

        for (size_t i = 0; i < 3; i++) {
            int i1 = point_id[i];
            for (int j = 0; j < 3; j++) {
                int j1 = point_id[j];
                global_A[i1 * num_of_points + j1] += local_A[i * 3 + j];
                global_F[i1] += local_F[i * 3 + j] * F;
            }
        }
    }

    for (size_t i = 0; i < num_of_points; i++) {
        if (out.pointmarkerlist[i] == 1) {
            for (size_t j = 0; j < num_of_points; j++) {
                global_F[j] -= global_A[j * num_of_points + i] * boundary_condition(i, &out);
                global_A[i * num_of_points + j] = 0;
                global_A[j * num_of_points + i] = 0;
            }
            global_A[i * num_of_points + i] = 1.0;
            global_F[i] = boundary_condition(i, &out);
        }
    }


    /*************************************************************************
     ************************ Solve Ax = b ***********************************
     *************************************************************************/
    int ipiv[num_of_points];
    int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR,
                             num_of_points,
                             1,
                             global_A,
                             num_of_points,
                             ipiv,
                             global_F,
                             1);
    if (info > 0) {
        printf("The factorization has been completed, but the factor U is\n");
        printf("exactly singular, so the solution could not be computed.");
        // todo goto mem free
        return 1;
    } else if (info < 0) {
        printf("%d-th argument had an illegal value\n", -info);
        // todo goto mem free
        return 1;
    }

    printf("\nSolution in point vs exact solution\n");
    for (size_t i = 0; i < num_of_points; i++) {
        double temp = exact_solution(i, &out);
        printf(
            "%lf  \t  %lf  \t  %lf\n", temp, global_F[i], temp - global_F[i]);
    }

    /* free all Triangle data structures */
    free(out.edgelist);
    free(out.edgemarkerlist);
    free(out.neighborlist);
    free(out.pointlist);
    free(out.trianglelist);
    free(out.segmentlist);
    free(out.segmentmarkerlist);
    free(out.pointmarkerlist);

    /* free all input data structures */
    free(in.regionlist);
    free(in.pointlist);

    /* free all matrix arrays */
    free(global_F);
    free(local_F);

    return 0;
}

static void generate_mesh(struct triangulateio *in, struct triangulateio *out)
{
    in->numberofpoints = 4;
    in->numberofpointattributes = 0;
    in->pointattributelist = (REAL *) NULL;
    in->pointlist = malloc(in->numberofpoints * 2 * sizeof(REAL));

    in->pointlist[0] = 0.0;
    in->pointlist[1] = 0.0;

    in->pointlist[2] = 0.0;
    in->pointlist[3] = 1.0;

    in->pointlist[4] = 1.0;
    in->pointlist[5] = 1.0;

    in->pointlist[6] = 1.0;
    in->pointlist[7] = 0.0;

    /*
     * "If you want Triangle to determine for you which vertices and edges are
     *  on the boundary, assign them the boundary marker zero (or use no
     *   markers at all) in your input files. In the output files, all boundary
     *   vertices, edges, and segments will be assigned the value one."
     */
    in->pointmarkerlist = (int *) NULL;

    in->numberofsegments = 0;
    in->numberofholes = 0;
    in->numberofregions = 1;
    in->regionlist = malloc(in->numberofregions * 4 * sizeof(REAL));
    in->regionlist[0] = 0.5;
    in->regionlist[1] = 0.5;
    in->regionlist[2] = 7.0; /* Regional attribute (for whole mesh). */
    in->regionlist[3] = 0.1; /* Area constraint. */

    out->pointlist = (REAL *) NULL;
    out->pointmarkerlist = (int *) NULL;
    out->trianglelist = (int *) NULL;
    out->triangleattributelist = (REAL *) NULL;
    out->neighborlist = (int *) NULL;
    out->segmentlist = (int *) NULL;
    out->segmentmarkerlist = (int *) NULL;
    out->edgelist = (int *) NULL;
    out->edgemarkerlist = (int *) NULL;

    /* Triangulate the points. read and write a PSLG (p), preserve the
     * convex hull (c), zero index (z), produce edge list (e), produce
     * neighbor list (n), be quiet (Q), area should be less than 0.1 (a.1)
     * and generate triangle for FEM (q).
     */
    triangulate("pzceQna.0001q", in, out, (struct triangulateio *) NULL);
}

static REAL exact_solution(uint point_id, struct triangulateio *out)
{
    REAL x = out->pointlist[2 * point_id];
    REAL y = out->pointlist[2 * point_id + 1];
    return 1.0 + x * x + 2.0 * y * y;
}

static REAL boundary_condition(uint point_id, struct triangulateio *out)
{
    REAL x = out->pointlist[2 * point_id];
    REAL y = out->pointlist[2 * point_id + 1];
    return 1.0 + x * x + 2.0 * y * y;
}
