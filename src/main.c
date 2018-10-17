#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

#include "triangle/triangle.h"

static void compute_local(int *point_id, REAL *pointlist, REAL (*matrix)[3]);
static void print_local_matrix(REAL *matrix);
static void print_global_matrix(REAL *matrix, int num_elements, bool num);
static void generate_mesh(struct triangulateio *io, struct triangulateio *out);

int main(void)
{
    struct triangulateio in, out;
    REAL(*l_matrix)[3] = malloc(3 * sizeof(*l_matrix));

    /* generate geometry */
    generate_mesh(&in, &out);

    REAL(*g_matrix)
    [out.numberofpoints] = malloc(out.numberofpoints * sizeof(*g_matrix));
    memset(g_matrix, 0, out.numberofpoints * sizeof(*g_matrix));

    for (int n = 0; n < out.numberoftriangles; n++) {
        // three points that describe triangle
        int point_id[3];
        memcpy(point_id, &out.trianglelist[n * 3], 3 * sizeof(int));
        compute_local(point_id, out.pointlist, l_matrix);

        for (int i = 0; i < 3; i++) {
            int i1 = point_id[i];
            for (int j = 0; j < 3; j++) {
                int j1 = point_id[j];
                g_matrix[i1][j1] += l_matrix[i][j];
            }
        }
    }
    print_global_matrix((REAL *) g_matrix, out.numberofpoints, false);

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
    free(g_matrix);
    free(l_matrix);

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

    /* Triangulate the points. Read and write a PSLG (p), preserve the
     * convex hull (c), zero index (z), produce edge list (e), produce
     * neighbor list (n), be quiet (Q), area should be less than 0.1 (a.1)
     * and generate triangle for FEM (q).
     */
    triangulate("pzceQna.1q", in, out, (struct triangulateio *) NULL);
}

static void compute_local(int *point_id, REAL *pointlist, REAL (*matrix)[3])
{
    REAL x1 = pointlist[2 * point_id[0]];
    REAL y1 = pointlist[2 * point_id[0] + 1];
    REAL x2 = pointlist[2 * point_id[1]];
    REAL y2 = pointlist[2 * point_id[1] + 1];
    REAL x3 = pointlist[2 * point_id[2]];
    REAL y3 = pointlist[2 * point_id[2] + 1];

    REAL dx23 = x2 - x3;
    REAL dy23 = y2 - y3;
    REAL dx31 = x3 - x1;
    REAL dy31 = y3 - y1;
    REAL dx12 = x1 - x2;
    REAL dy12 = y1 - y2;

    REAL area = 0.5 * (dx31 * dy12 - dy31 * dx12);

    matrix[0][0] = 0.25 * (dx23 * dx23 + dy23 * dy23) / area;
    matrix[0][1] = 0.25 * (dx23 * dx31 + dy23 * dy31) / area;
    matrix[0][2] = 0.25 * (dx23 * dx12 + dy23 * dy12) / area;

    matrix[1][0] = 0.25 * (dx31 * dx23 + dy31 * dy23) / area;
    matrix[1][1] = 0.25 * (dx31 * dx31 + dy31 * dy31) / area;
    matrix[1][2] = 0.25 * (dx31 * dx12 + dy31 * dy12) / area;

    matrix[2][0] = 0.25 * (dx12 * dx23 + dy12 * dy23) / area;
    matrix[2][1] = 0.25 * (dx12 * dx31 + dy12 * dy31) / area;
    matrix[2][2] = 0.25 * (dx12 * dx12 + dy12 * dy12) / area;
}

static void print_local_matrix(REAL *matrix)
{
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%lf ", matrix[i * 3 + j]);
        }
        printf("\n");
    }
}

static void print_global_matrix(REAL *matrix, int num_elements, bool num)
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
