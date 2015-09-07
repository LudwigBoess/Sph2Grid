/*Here we initialise the common variables*/
#ifndef COMMON_H

#include "allvars.h"
#include "schemes.h"

/* 
 * Global Variables
 */
struct InputParameters Param;
struct Snap_Properties Snap;
struct Particle_Data *P;
struct Gas_Data *Gas;
struct Grid_Data *Grid;

struct ParallelInfos ThisTask = { 0, 0, 0, {0}, 0, 0 };

#ifdef OPENMP
struct OpenMP_infos Omp = { 0, 0 };
#endif

/* Convert 3D to 1D index to access Grid */
inline size_t Idx(ptrdiff_t x, ptrdiff_t y, ptrdiff_t z)
{
    const int nGrid = Param.NGrid;

#ifdef PERIODIC			// This handles out of box indizes 
	while (x < 0)
		x += nGrid;

	while (x >= nGrid)
		x -= nGrid;

	while (y < 0)
		y += nGrid;
	
    while (y >= nGrid)
		y -= nGrid;

	while (z < 0)
		z += nGrid;
	
    while (z >= nGrid)
		z -= nGrid;
#endif // PERIODIC

	return (x * nGrid * nGrid + y * nGrid + z);
}

/* catch malloc errors */
inline void *my_malloc(size_t size)
{
	void *result = malloc(size);

	my_assert(result != NULL, "Could not allocate memory");

	return (result);
}

/* catch realloc errors */
inline void *my_realloc(void *ptr, size_t size)
{
	void *result = realloc(ptr, size);

	my_assert((result != NULL) || (size == 0), "Could not reallocate memory");

	return (result);
}

/* don't free a NULL ptr */
void my_free(void *ptr)
{
	if (!(ptr == NULL))
		free(ptr);

	return;
}

/* 
 * Reallocates the Particle structures
 * P and Gas. Takes the relative change
 * as argument, not the total number.
 * Add or Remove via sign argument. 
 */
extern void
Reallocate_P(long long partTotal, long long nPart[N_part_types], int sign)
{
	int type = 0;

	if (!partTotal)
		return;

	ThisTask.PartTotal += sign * llabs(partTotal);
	for (type = 0; type < N_part_types; type++) {
		ThisTask.Npart[type] += sign * llabs(nPart[type]);

		my_assert(ThisTask.Npart[type] >= 0,
			  "Can't allocate negative particles");
	}

	P = my_realloc((void *) P, sizeof(*P) * ThisTask.PartTotal);

	Gas = my_realloc((void *) Gas, sizeof(*Gas) * ThisTask.Npart[0]);

	if (ThisTask.PartTotal == 0)
		printf("No particles on processor <%d>\n", ThisTask.Rank);

	return;
}

/* 
 * Parallel data distribution along x-axis
 * For weird NTask nGrid constellations the 
 * fftw3 decomposition throws tasks out. As we
 * prefer memory over speed we always include
 * all tasks. See fftw3 docu for fftw_mpi_local_size_many 
 */
void distribute_grid()
{
	const double nGrid = (double) Param.NGrid;

	if (!ThisTask.Rank)
		printf("Parallel grid decomposition ... ");

	ThisTask.Nx = floor(nGrid / ThisTask.NTask);

	ThisTask.X_start = ThisTask.Rank * ThisTask.Nx;

	if (ThisTask.Rank == ThisTask.NTask - 1)
		ThisTask.Nx += nGrid - ThisTask.Nx * ThisTask.NTask;

	if (Param.NGrid % ThisTask.NTask && !ThisTask.Rank)
		printf("\n   Distribution not ideal !\n");

	ThisTask.NDouble = nGrid * nGrid * ThisTask.Nx;

	size_t isrc = Idx(ThisTask.X_start, 0, 0);

	size_t nBytes = ThisTask.NDouble * sizeof(*Grid);

	memmove(&(Grid[0]), &(Grid[isrc]), nBytes);

	Grid = my_realloc(Grid, nBytes); // now the Grid is stored distributed

	MPI_Barrier(MPI_COMM_WORLD);

	if (!ThisTask.Rank)
		printf("done \n\n");

	return;
}

/* Error Handling */
void my_assert_full(int expr, const char *errmsg,
		    const char *func, const char *file, int line)
{
	if (expr)
		return;

	fprintf(stderr, "\n\n(%d) %s:%d : %s : %s\n\n\n",
		ThisTask.Rank, file, line, func, errmsg);

	/* end program */
	MPI_Abort(MPI_COMM_WORLD, -1);
	exit(-1);
}

void null_final_operations()
{
	return;
}

/* Normalisation of the grid */
void divide_npart_final_operations()
{				
	enum gridfields block; // average over particles 
	size_t i, j;

	const size_t nDouble = sizeof(*Grid) / sizeof(double);

	for (block = (enum gridfields) 0; block < GRID_LASTENTRY; block++) {

		set_gridfield_prop(block);

		if (!GFld.rank)		// == not present 
			continue;

		if (block == GRID_MASS)	// in SPH mass is additive 
			continue;

		if (block == GRID_NPART) // don't be silly 
			continue;

		for (i = 0; i < NGrid3; i++) {

			if (!(long) (Grid[i].Npart))
				continue;

			for (j = 0; j < GFld.rank; j++)
				*(GFld.ptr + i * nDouble + j) /= Grid[i].Npart;
		}
	}

	return;
}

#endif
