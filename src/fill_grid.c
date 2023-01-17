#include "allvars.h"
#include "schemes.h"

#define FILTER 6		/* Filter size in cells -1  */

static inline void limit_idx_to_box(struct Cell_boundary_indices *);
static inline void distribute_particle(const float, const size_t, 
        const size_t);

static void reduce_grid();
static void compute_local_turbulence();

void fill_grid()
{
	start_timing(CPU_FILL);

	const float halfGridSize = Param.GridSize * 0.5;

	if (!ThisTask.Rank) 
		printf("Filling  Grid ... \n");

#pragma omp parallel for // in cell units
	for (int ipart = 0; ipart < ThisTask.PartTotal; ipart++) { 
		
		float px = (P[ipart].Pos[0] + halfGridSize) * Kpc2cell;
		float py = (P[ipart].Pos[1] + halfGridSize) * Kpc2cell;
		float pz = (P[ipart].Pos[2] + halfGridSize) * Kpc2cell;

		float psize = WMODE_PART_SIZE(ipart) * Kpc2cell;

	    struct Cell_boundary_indices idx = { 0 }; // is send around a lot 

		WMODE_EDGE_IDX(px, py, pz, psize, &idx); // find particle edges 
	
		limit_idx_to_box(&idx);

		/* distribute, idx < 0 possible */
		for (ptrdiff_t xidx = idx.Xmin; xidx < idx.Xmax; xidx++) {

			float x = xidx + 0.5;

			for (ptrdiff_t yidx = idx.Ymin; yidx < idx.Ymax; yidx++) {

				float y = yidx + 0.5;

				for (ptrdiff_t zidx = idx.Zmin; zidx < idx.Zmax; zidx++) {

					float z = zidx + 0.5;

					size_t gridIdx = Idx(xidx, yidx, zidx);	//handles PERIODIC
					
                    float weight = WMODE_WEIGHTS(ipart, x, y, z, px, py, pz);

					if (!weight)
						continue;
					
                    distribute_particle(weight, ipart, gridIdx);
				}
			}
		}
	}
	
	stop_timing(CPU_FILL);

	reduce_grid();

	WMODE_FINAL_OPS();

	compute_local_turbulence();

	return;
}

/* this one is optimised out */
static inline void distribute_particle(const float weight, const size_t ipart, 
        const size_t idx)
{
    Grid[idx].Mass += P[ipart].Mass * weight; 
    
    Grid[idx].Rho += Gas[ipart].Rho * weight;

#ifdef VEL
	Grid[idx].Vel[0] += P[ipart].Vel[0] * weight;
	Grid[idx].Vel[1] += P[ipart].Vel[1] * weight;
	Grid[idx].Vel[2] += P[ipart].Vel[2] * weight;
#endif /* VEL */

#ifdef MOMENTUM
	Grid[idx].Mom[0] += P[ipart].Vel[0] * P[ipart].Mass * weight;
	Grid[idx].Mom[1] += P[ipart].Vel[1] * P[ipart].Mass * weight;
	Grid[idx].Mom[2] += P[ipart].Vel[2] * P[ipart].Mass * weight;
#endif /* MOMENTUM */

#ifdef INTENERGY
	Grid[idx].U += Gas[ipart].U * weight;
#endif /* INTENERGY */

#ifdef BFLD
	Grid[idx].Bfld[0] += Gas[ipart].Bfld[0] * weight;
	Grid[idx].Bfld[1] += Gas[ipart].Bfld[1] * weight;
	Grid[idx].Bfld[2] += Gas[ipart].Bfld[2] * weight;
#endif // BFLD

#ifdef SCALAR_BFLD
    Grid[idx].B += length3(Gas[ipart].Bfld) * weight;
#endif // SCALAR_BFLD

#ifdef DENSVEL			 
	double tmp = pow(Gas[ipart].Rho, 1. / 3.);
	
    Grid[idx].DensVel[0] += tmp * P[ipart].Vel[0] * weight;
	Grid[idx].DensVel[1] += tmp * P[ipart].Vel[1] * weight;
	Grid[idx].DensVel[2] += tmp * P[ipart].Vel[2] * weight;
#endif /* DENSVEL */

#ifdef SCALARVEL
	Grid[idx].VelScalar += length3(P[ipart].Vel) * weight;
#endif /* SCALARVEL */

#ifdef VTURB
	Grid[idx].VTurb += Gas[ipart].VTurb * weight;
#endif

#ifdef MACH
	Grid[idx].Mach += Gas[ipart].Mach * weight;
#endif
	
    Grid[idx].Npart++;

    return;
}

/* ensure  0 <= idx < Ngrid  */
static inline void limit_idx_to_box(struct Cell_boundary_indices *idx)
{
#ifndef PERIODIC
	idx->Xmin = fmax(idx->Xmin, 0);
	idx->Xmax = fmin(idx->Xmax, Param.NGrid-1);

	idx->Ymin = fmax(idx->Ymin, 0);
	idx->Ymax = fmin(idx->Ymax, Param.NGrid-1);

	idx->Zmin = fmax(idx->Zmin, 0);
	idx->Zmax = fmin(idx->Zmax, Param.NGrid-1);
#endif // with PERIODIC  out of box idx are handled in Idx(i,j,k) 

	return;
}

/* Add up all Quantities over all CPUs */
static void reduce_grid()
{
    if (ThisTask.NTask == 1)
        return;

    start_timing(CPU_REDUCE);

	double * restrict sendRecvBuf = my_malloc(NGrid3 * sizeof(*sendRecvBuf));

	const size_t nDouble = sizeof(*Grid) / sizeof(double);

	for (int block = 0; block < GRID_LASTENTRY; block++) { // Reduce and write 

		set_gridfield_prop((enum gridfields) block);

		if (!GFld.rank)
			continue;

		for (size_t j = 0; j < GFld.rank; j++) {	// component-wise 

			for (size_t i = 0; i < NGrid3; i++) 
				sendRecvBuf[i] = *(GFld.ptr + i*nDouble + j);

			MPI_Allreduce(MPI_IN_PLACE, sendRecvBuf, NGrid3,
				      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			for (size_t i = 0; i < NGrid3; i++)
				*(GFld.ptr + i*nDouble + j) = sendRecvBuf[i];
		}
	}

	my_free(sendRecvBuf);

	stop_timing(CPU_REDUCE);

	return;
}

/* This is a second particle loop to compute
 * derived quantities */
static void compute_local_turbulence()
{
#if defined(VELDISPERSION) || defined(VELFILTERED)
	size_t ipart;
	double px, py, pz;
	ptrdiff_t xidx, yidx, zidx;
	size_t gridIdx, sumIdx = 0;

	struct Cell_boundary_indices idx;

	if (!ThisTask.Rank)
		printf("Makeing local velocity dispersion \n\n");

	double halfGridSize = Param.GridSize * 0.5;


#ifdef VELFILTERED
	if (!ThisTask.Rank)
		printf("Filter size is %d cells, or %g kpc \n", FILTER,
		       FILTER * Param.GridSize / Param.NGrid);
#endif
#ifdef OPENMP
#pragma omp parallel for \
        shared(Grid) \
        private(ipart, px,py,pz,xidx,yidx,zidx, gridIdx, idx)
#endif
	for (ipart = 0; ipart < ThisTask.PartTotal; ipart++) {
		/* convert to cell units */
		px = (P[ipart].Pos[0] + halfGridSize) * Kpc2cell;
		py = (P[ipart].Pos[1] + halfGridSize) * Kpc2cell;
		pz = (P[ipart].Pos[2] + halfGridSize) * Kpc2cell;

#ifdef VELDISPERSION
		/* NGP */
		idx.Xmin = floor(px);
		idx.Ymin = floor(py);
		idx.Ymax = idx.Ymin + 1;
		idx.Zmin = floor(pz);
		idx.Zmax = idx.Zmin + 1;
#endif
#ifdef VELFILTERED
		idx.Xmin = floor(px) - FILTER / 2;
		idx.Xmax = idx.Xmin + FILTER;

		idx.Ymin = floor(py) - FILTER / 2;
		idx.Ymax = idx.Ymin + FILTER;
		
        idx.Zmin = floor(pz) - FILTER / 2;
		idx.Zmax = idx.Zmin + FILTER;
		
        sumIdx = gridIdx = Idx(floor(px), floor(py), floor(pz));
#endif
		for (xidx = idx.Xmin; xidx < idx.Xmax; xidx++) {
			for (yidx = idx.Ymin; yidx < idx.Ymax; yidx++) {
				for (zidx = idx.Zmin; zidx < idx.Zmax; zidx++) {

					gridIdx = Idx(xidx, yidx, zidx);	/* handles PERIODIC */

#ifdef VELDISPERSION
					Grid[gridIdx].Sigma[0] +=
					    p2(P[ipart].Vel[0] - Grid[gridIdx].Vel[0]);
					Grid[gridIdx].Sigma[1] +=
					    p2(P[ipart].Vel[1] - Grid[gridIdx].Vel[1]);
					Grid[gridIdx].Sigma[2] +=
					    p2(P[ipart].Vel[2] - Grid[gridIdx].Vel[2]);
#endif
#ifdef VELFILTERED		/* this becomes the mean */
					Grid[sumIdx].VelFiltered[0] +=
					    Grid[gridIdx].Vel[0] / p3(FILTER + 1);
					Grid[sumIdx].VelFiltered[1] +=
					    Grid[gridIdx].Vel[1] / p3(FILTER + 1);
					Grid[sumIdx].VelFiltered[2] +=
					    Grid[gridIdx].Vel[2] / p3(FILTER + 1);
#endif
				}
			}
		}
	}

	for (gridIdx = 0; gridIdx < NGrid3; gridIdx++) {
#ifdef VELFILTERED
		Grid[gridIdx].VelFiltered[0] = Grid[gridIdx].Vel[0] 
			- Grid[gridIdx].VelFiltered[0];
		Grid[gridIdx].VelFiltered[1] = Grid[gridIdx].Vel[1] 
			- Grid[gridIdx].VelFiltered[1];
		Grid[gridIdx].VelFiltered[2] = Grid[gridIdx].Vel[2] 
			- Grid[gridIdx].VelFiltered[2];
#endif
#ifdef VELDISPERSION
		if (!Grid[gridIdx].Npart)
			continue;

		Grid[gridIdx].Sigma[0] = sqrt(Grid[gridIdx].Sigma[0] 
				/ Grid[gridIdx].Npart);
		Grid[gridIdx].Sigma[1] = sqrt(Grid[gridIdx].Sigma[1] 
				/ Grid[gridIdx].Npart);
		Grid[gridIdx].Sigma[2] = sqrt(Grid[gridIdx].Sigma[2] 
				/ Grid[gridIdx].Npart);
#endif

	}

#endif // defined(VELDISPERSION) || defined(VELFILTERED 
	return;
}
