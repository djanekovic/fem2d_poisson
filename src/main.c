#include <stdio.h>
#include <stdlib.h>

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

#include "triangle/triangle.h"

static void report(io, markers, reporttriangles, reportneighbors, reportsegments,
            reportedges, reportnorms)
struct triangulateio *io;
int markers;
int reporttriangles;
int reportneighbors;
int reportsegments;
int reportedges;
int reportnorms;
{
  int i, j;

  for (i = 0; i < io->numberofpoints; i++) {
    printf("Point %4d:", i);
    for (j = 0; j < 2; j++) {
      printf("  %.6g", io->pointlist[i * 2 + j]);
    }
    if (io->numberofpointattributes > 0) {
      printf("   attributes");
    }
    for (j = 0; j < io->numberofpointattributes; j++) {
      printf("  %.6g",
             io->pointattributelist[i * io->numberofpointattributes + j]);
    }
    if (markers) {
      printf("   marker %d\n", io->pointmarkerlist[i]);
    } else {
      printf("\n");
    }
  }
  printf("\n");

  if (reporttriangles || reportneighbors) {
    for (i = 0; i < io->numberoftriangles; i++) {
      if (reporttriangles) {
        printf("Triangle %4d points:", i);
        for (j = 0; j < io->numberofcorners; j++) {
          printf("  %4d", io->trianglelist[i * io->numberofcorners + j]);
        }
        if (io->numberoftriangleattributes > 0) {
          printf("   attributes");
        }
        for (j = 0; j < io->numberoftriangleattributes; j++) {
          printf("  %.6g", io->triangleattributelist[i *
                                         io->numberoftriangleattributes + j]);
        }
        printf("\n");
      }
      if (reportneighbors) {
        printf("Triangle %4d neighbors:", i);
        for (j = 0; j < 3; j++) {
          printf("  %4d", io->neighborlist[i * 3 + j]);
        }
        printf("\n");
      }
    }
    printf("\n");
  }

  if (reportsegments) {
    for (i = 0; i < io->numberofsegments; i++) {
      printf("Segment %4d points:", i);
      for (j = 0; j < 2; j++) {
        printf("  %4d", io->segmentlist[i * 2 + j]);
      }
      if (markers) {
        printf("   marker %d\n", io->segmentmarkerlist[i]);
      } else {
        printf("\n");
      }
    }
    printf("\n");
  }

  if (reportedges) {
    for (i = 0; i < io->numberofedges; i++) {
      printf("Edge %4d points:", i);
      for (j = 0; j < 2; j++) {
        printf("  %4d", io->edgelist[i * 2 + j]);
      }
      if (reportnorms && (io->edgelist[i * 2 + 1] == -1)) {
        for (j = 0; j < 2; j++) {
          printf("  %.6g", io->normlist[i * 2 + j]);
        }
      }
      if (markers) {
        printf("   marker %d\n", io->edgemarkerlist[i]);
      } else {
        printf("\n");
      }
    }
    printf("\n");
  }
}

int main(void) {

    /* generate geometry */
    struct triangulateio in, out;

    in.numberofpoints = 4;
    in.numberofpointattributes = 0;
    in.pointattributelist = (REAL *) NULL;
    in.pointlist = malloc(in.numberofpoints * 2 * sizeof(REAL));

    in.pointlist[0] = 0.0;
    in.pointlist[1] = 0.0;

    in.pointlist[2] = 0.0;
    in.pointlist[3] = 1.0;

    in.pointlist[4] = 1.0;
    in.pointlist[5] = 1.0;

    in.pointlist[6] = 1.0;
    in.pointlist[7] = 0.0;

    /*
     * "If you want Triangle to determine for you which vertices and edges are
     *  on the boundary, assign them the boundary marker zero (or use no
     *   markers at all) in your input files. In the output files, all boundary
     *   vertices, edges, and segments will be assigned the value one."
     */
    in.pointmarkerlist = (int *) NULL;

    in.numberofsegments = 0;
    in.numberofholes = 0;
    in.numberofregions = 1;
    in.regionlist = malloc(in.numberofregions * 4 * sizeof(REAL));
    in.regionlist[0] = 0.5;
    in.regionlist[1] = 0.5;
    in.regionlist[2] = 7.0;            /* Regional attribute (for whole mesh). */
    in.regionlist[3] = 0.1;          /* Area constraint. */

    out.pointlist = (REAL *) NULL;
    out.pointmarkerlist = (int *) NULL;
    out.trianglelist = (int *) NULL;
    out.triangleattributelist = (REAL *) NULL;
    out.neighborlist = (int *) NULL;
    out.segmentlist = (int *) NULL;
    out.segmentmarkerlist = (int *) NULL;
    out.edgelist = (int *) NULL;
    out.edgemarkerlist = (int *) NULL;

    /* Triangulate the points. Read and write a PSLG (p), preserve the
     * convex hull (c), zero index (z), produce edge list (e), produce
     * neighbor list (n), be quiet (Q), area should be less than 0.1 (a.1)
     * and generate triangle for FEM (q).
     */
    triangulate("pzceQna.1q", &in, &out, (struct triangulateio *) NULL);

    //report(&out, 1, 1, 1, 1, 1, 0);

    return 0;
}
