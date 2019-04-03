/**
 * Finite element method solver. This solver solves laplace(u) = 6 on square
 * [0, 1]x[0, 1]. Boundary conditions are set to be 1 + x^2 + 2y^2.
 *
 * Exact solution of this problem is 1 + x^2 + 2y^2. We will test against this
 * solution in the end.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <lapacke.h>
#include <triangle.h>

#define F -6.0

#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

static void compute_local_a(uint *point_id, REAL *pointlist, REAL *matrix);
static void compute_local_f(uint *point_id, REAL *pointlist, REAL *matrix);
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
        compute_local_a(point_id, out.pointlist, local_A);
        compute_local_f(point_id, out.pointlist, local_F);

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
                global_F[j] -= global_A[j * num_of_points + i]
                    * boundary_condition(i, &out);
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
        goto mem_free;
    } else if (info < 0) {
        printf("%d-th argument had an illegal value\n", -info);
        goto mem_free;
    }

    //TODO: implement max error
    printf("\nSolution in point vs exact solution\n");
    for (size_t i = 0; i < num_of_points; i++) {
        double temp = exact_solution(i, &out);
        printf(
            "%lf  \t  %lf  \t  %lf\n", temp, global_F[i], temp - global_F[i]);
    }

mem_free:
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
    free(global_A);
    free(local_A);
    free(global_F);
    free(local_F);

    return info;
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
    triangulate("pzceQna.1q", in, out, (struct triangulateio *) NULL);
}

static void compute_local_a(uint *point_id, REAL *pointlist, REAL *matrix)
{
        REAL dx23 = pointlist[2 * point_id[1]] - pointlist[2 * point_id[2]];
        REAL dx31 = pointlist[2 * point_id[2]] - pointlist[2 * point_id[0]];
        REAL dx12 = pointlist[2 * point_id[0]] - pointlist[2 * point_id[1]];
        REAL dy23 = pointlist[2 * point_id[1] + 1] - pointlist[2 * point_id[2] + 1];
        REAL dy31 = pointlist[2 * point_id[2] + 1] - pointlist[2 * point_id[0] + 1];
        REAL dy12 = pointlist[2 * point_id[0] + 1] - pointlist[2 * point_id[1] + 1];;

        REAL area = 0.5 * (dx31 * dy12 - dy31 * dx12);
        REAL _tmp_mult = 0.25 / area;
        matrix[0] = _tmp_mult * (dx23 * dx23 + dy23 * dy23);
        matrix[1] = _tmp_mult * (dx23 * dx31 + dy23 * dy31);
        matrix[2] = _tmp_mult * (dx23 * dx12 + dy23 * dy12);
        matrix[3] = _tmp_mult * (dx31 * dx23 + dy31 * dy23);
        matrix[4] = _tmp_mult * (dx31 * dx31 + dy31 * dy31);
        matrix[5] = _tmp_mult * (dx31 * dx12 + dy31 * dy12);
        matrix[6] = _tmp_mult * (dx12 * dx23 + dy12 * dy23);
        matrix[7] = _tmp_mult * (dx12 * dx31 + dy12 * dy31);
        matrix[8] = _tmp_mult * (dx12 * dx12 + dy12 * dy12);
}

static inline void compute_local_f(uint *point_id, REAL *pointlist, REAL *matrix)
{
        REAL dx31 = pointlist[2 * point_id[2]] - pointlist[2 * point_id[0]];
        REAL dx12 = pointlist[2 * point_id[0]] - pointlist[2 * point_id[1]];
        REAL dy31 = pointlist[2 * point_id[2] + 1] - pointlist[2 * point_id[0] + 1];
        REAL dy12 = pointlist[2 * point_id[0] + 1] - pointlist[2 * point_id[1] + 1];

        REAL area = 0.5 * (dx31 * dy12 - dy31 * dx12);
        REAL c_diag = area / 6.0;
        REAL c_off = area / 12.0;

        matrix[0] = c_diag; matrix[1] = c_off;  matrix[2] = c_off;
        matrix[3] = c_off;  matrix[4] = c_diag; matrix[5] = c_off;
        matrix[6] = c_off;  matrix[7] = c_off;  matrix[8] = c_diag;
}

static inline REAL exact_solution(uint point_id, struct triangulateio *out)
{
    REAL x = out->pointlist[2 * point_id];
    REAL y = out->pointlist[2 * point_id + 1];
    return (1.0 + x * x + 2.0 * y * y);
}

static inline REAL boundary_condition(uint point_id, struct triangulateio *out)
{
    REAL x = out->pointlist[2 * point_id];
    REAL y = out->pointlist[2 * point_id + 1];
    return (1.0 + x * x + 2.0 * y * y);
}
