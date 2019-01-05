#include <petscdm.h>
#include <petscdraw.h>
#include <petscdmda.h>
#include <petscdmplex.h>

/**
 * Finite element method solver. This solver solves laplace(u) = 6 on square
 * [0, 1]x[0, 1]. Boundary conditions are set to be 1 + x^2 + 2y^2.
 *
 * Exact solution of this problem is 1 + x^2 + 2y^2. We will test against this
 * solution in the end.
 *
 * Weak formulation is: inner(grad(u), grad(v))*dS = -6 *v*dS
 */

static PetscInt size, rank;
static const char help[] = "Solving poisson equation in 2D\n\n";

struct options {
    PetscInt dim;           /* mesh dimensionality, defaults to 2            */
    PetscInt faces[3];      /* Number of faces in xyz, default (2, 2, 0)     */

    PetscInt c_start;       /* Cell start/end                                */
    PetscInt c_end;
    PetscInt e_start;       /* Edge start/end                                */
    PetscInt e_end;
    PetscInt v_start;       /* Vertex start/end                              */
    PetscInt v_end;
};

static PetscErrorCode boundary_condition(const PetscReal x[], PetscScalar *u) {
    *u = x[0] * x[0] + 2 * x[1] * x[1] + 1;
    return 0;
}

static void set_default(struct options **ctx)
{
    (*ctx)->dim = 2;

    (*ctx)->faces[0] = 2;
    (*ctx)->faces[1] = 2;
    (*ctx)->faces[2] = 0;
}

static PetscErrorCode read_cli_options(struct options *ctx)
{
    PetscErrorCode ierr;

    set_default(&ctx);

    // Number of dimensions
    ierr = PetscOptionsGetInt(NULL, NULL, "-dim", &(ctx->dim), NULL); CHKERRQ(ierr);

    if ( ctx->dim == 3) {
        PetscPrintf(PETSC_COMM_WORLD, "3D not yet supported\n");
        PetscPrintf(PETSC_COMM_WORLD, "Exiting...\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "3D not yet implemented\n");
    }

    // Number of faces in x, y, z
    ierr = PetscOptionsGetInt(NULL, NULL, "-x", &(ctx->faces[0]), NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-y", &(ctx->faces[1]), NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-z", &(ctx->faces[2]), NULL); CHKERRQ(ierr);

    return 0;
}

static PetscErrorCode create_mesh(struct options ctx, DM *dm)
{
    PetscInt id=1;
    DMLabel label;
    PetscErrorCode ierr;

    /* Using world communicator, spatial dimension is dim, generate simplices, */
    ierr = DMPlexCreateBoxMesh(PETSC_COMM_WORLD, ctx.dim, PETSC_TRUE, NULL, NULL, NULL, NULL, PETSC_TRUE, dm);
    CHKERRQ(ierr);

    DMCreateLabel(*dm, "boundary");
    ierr = DMGetLabel(*dm, "boundary", &label); CHKERRQ(ierr);
    ierr = DMPlexMarkBoundaryFaces(*dm, 1, label); CHKERRQ(ierr);
    ierr = DMPlexLabelComplete(*dm, label); CHKERRQ(ierr);
    DMAddBoundary(*dm, DM_BC_ESSENTIAL, "wall", "boundary", 0, 0, NULL, (void (*)()) boundary_condition, 1, &id, NULL);

    return 0;
}

static PetscErrorCode get_mesh_attributes(struct options *ctx, DM dm)
{
    PetscErrorCode ierr;
    PetscInt start, end;

    ierr = DMPlexGetHeightStratum(dm, 0, &start, &end);     /*  cells        */
    CHKERRQ(ierr);
    ctx->c_start = start;
    ctx->c_end = end;

    ierr = DMPlexGetHeightStratum(dm, 1, &start, &end);     /*  edges        */
    CHKERRQ(ierr);
    ctx->e_start = start;
    ctx->e_end = end;

    ierr = DMPlexGetHeightStratum(dm, 2, &start, &end);     /*  vertices     */
    CHKERRQ(ierr);
    ctx->v_start = start;
    ctx->v_end = end;

    return 0;
}
int main (int argc, char **argv)
{
    struct options ctx;
    Vec coordinates;
    DM dm;

    PetscErrorCode ierr = PetscInitialize(&argc, &argv, (char *) NULL, help);
    if (ierr) {
        return ierr;
    }

    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = read_cli_options(&ctx);
    ierr = create_mesh(ctx, &dm);
    ierr = get_mesh_attributes(&ctx, dm);

    DMGetCoordinates(dm, &coordinates);
    VecView(coordinates, PETSC_VIEWER_STDOUT_WORLD);

    PetscFinalize();
    return 0;
}
