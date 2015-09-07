#include "allvars.h"

#ifdef FOURIERTRANSFORM

static void clean_FFT_dataspace(double *, complex double *);
static void prepare_FFT_input(int, double *);
static void set_FFT_normalisation(complex double *);
static void move_FFT_output(int, complex double *);

static ptrdiff_t nx, ny, nz, nz_padded;

/* 
 * Compute FFT of Grid in parallel -
 * the Hermitian symmetry is used to store only 
 * the non-redundant data. The x-dimension is
 * split over all tasks. 
 */
void fft_grid()
{
	int block, comp;
	ptrdiff_t fftw_alloc_size, n[3] = { 0 };
	double *rdata_in;
	complex double *cdata_out;	/* C99 */
	fftw_plan plan;

	start_timing(CPU_FFT);

	if (!ThisTask.Rank)
		printf("Starting parallel FFT ... \n");

	fftw_mpi_init();

#if defined(OPENMP) && defined(OMP_FFTW3)
	fftw_init_threads();
	fftw_plan_with_nthreads(Omp.NThreads);
#endif

	/* number of cells */
	nx = ny = nz = Param.NGrid;
	nz_padded = 2 * (nz / 2 + 1);	/* FFTW3 Manual, sect. 2.4 */

	/* always include all CPUs. see distribute_grid() */
	n[0] = nx;
	n[1] = ny;
	n[2] = nz_padded;

	fftw_alloc_size = fftw_mpi_local_size_many(3, n, 1, 
			floor(Param.NGrid / ThisTask.NTask), MPI_COMM_WORLD,
			&ThisTask.FFT_nx, &ThisTask.FFT_x_start);

#ifdef HACKFFTW3		/* Ugly */
	// ThisTask.FFT_nx = ThisTask.Nx;
	printf("HACKFFT Rank %d NX %zu x_start %02zu n[0] %zd n[1] %zd n[2] %zd \n",
	       ThisTask.Rank, ThisTask.FFT_nx, ThisTask.FFT_x_start, n[0], n[1], n[2]);
#endif

	my_assert(ThisTask.Nx == ThisTask.FFT_nx
		  && ThisTask.X_start == ThisTask.FFT_x_start,
		  "FFT data distribution wrong");

	/* allocate work arrays - sect. 6.5 of FFTW3 Manual */
	rdata_in = fftw_alloc_real(2 * fftw_alloc_size);
	cdata_out = fftw_alloc_complex(fftw_alloc_size);

	/* make FFTW plan */
	plan = fftw_mpi_plan_dft_r2c_3d(nx, ny, nz, rdata_in, cdata_out,
					MPI_COMM_WORLD, FFTW_ESTIMATE);

	/* Actual size of local FFT data */
	ThisTask.FFT_nComplex = ThisTask.FFT_nx * ny * nz_padded / 2;

	for (block = 0; block < GRID_LASTENTRY; block++) {

		set_gridfield_prop((enum gridfields) block);	/* set global GFld var */

		if (!GFld.rank || !GFld.FTptr)	/* not present */
			continue;

		if (!ThisTask.Rank)
			printf("    %s\n", GFld.name);

		/* loop over components */
		for (comp = 0; comp < GFld.rank; comp++) {
			clean_FFT_dataspace(rdata_in, cdata_out);

			prepare_FFT_input(comp, rdata_in);

			fftw_execute(plan);

			set_FFT_normalisation(cdata_out);

			move_FFT_output(comp, cdata_out);
		}
	}
	if (!ThisTask.Rank)
		printf("\n");

/*
    long i,j,k,ii,jj,kk, idxc,idx;
    double  realo, imo, kvec;
    
    make_kGrid();

    printf("\n");

    if (ThisTask.NTask == 1)
    for (i=0; i<ThisTask.FFT_nx; i++)
    for (j=0; j<ny; j++)
    for (k=0; k<nz; k++){

        ii = i;
        jj = j;
        kk = k;

        if (k > nz/2){
            kk = nz - k;
            
            if ( j != 0 )
                jj = ny - j;

            if ( i != 0 )
                ii = nx - i;

        }

        idxc = (ii*ny + jj) * (nz_padded/2) + kk;

        realo = creal((FTGrid.Vel[0])[idxc]);
        imo = cimag((FTGrid.Vel[0])[idxc]);

        if ( k > nz/2)
            imo *= -1;

        kvec = sqrt( pow((FTGrid.kVector[0])[idxc],2) 
                + pow((FTGrid.kVector[1])[idxc],2)
                + pow((FTGrid.kVector[2])[idxc],2)  );

        printf("%d %d %d | %d %d %d | %d | %g %g %g \n",
                i,j,k, ii,jj,kk, idxc, realo, imo, kvec);

    }
MPI_Barrier(MPI_COMM_WORLD);
*/

	fftw_destroy_plan(plan);
	fftw_mpi_cleanup();
	fftw_free(rdata_in);
	fftw_free(cdata_out);
#if defined(OPENMP) && defined(OMP_FFTW3)
	fftw_cleanup_threads();
#endif

	stop_timing(CPU_FFT);

	return;

}

/* init/reset FFT work arrays */
void clean_FFT_dataspace(double *rdata_in, complex double *cdata_out)
{
	size_t i;

	for (i = 0; i < 2 * ThisTask.FFT_nComplex; i++)
		rdata_in[i] = 0;

	for (i = 0; i < ThisTask.FFT_nComplex; i++)
		cdata_out[i] = 0 + I * 0;

	return;
}

/* copy slab of grid to data needed by this task */
void prepare_FFT_input(int comp, double *data)
{
	ptrdiff_t i, j, k, grid_idx, data_idx, nDouble;

	nDouble = sizeof(*Grid) / sizeof(double);

	for (i = 0; i < ThisTask.FFT_nx; ++i) {
		for (j = 0; j < ny; ++j) {
			for (k = 0; k < nz; ++k) {
				grid_idx = (i * ny + j) * nz + k;
				data_idx = (i * ny + j) * nz_padded + k;

				data[data_idx] =
				    *(GFld.ptr + grid_idx * nDouble +
				      comp);
			}
		}
	}

	return;
}

/* Write FTGrid struct */
void move_FFT_output(int comp, complex double *data)
{
	size_t i, nBytes;

	nBytes = ThisTask.FFT_nComplex * sizeof(complex double);
	GFld.FTptr[comp] = my_malloc(nBytes);

	for (i = 0; i < ThisTask.FFT_nComplex; i++)
		(GFld.FTptr[comp])[i] = data[i];

	return;
}

/* There are just too many conventions here ... */
void set_FFT_normalisation(complex double *data)
{
	size_t i;
	double norm = 1;

#if   defined(FFT_IDL_NORM)
	norm = 1.0 / (nx * ny * nz);
#elif defined(FFT_SYM_NORM)
	norm = 1.0 / sqrt(nx * ny * nz);
#endif

	for (i = 0; i < ThisTask.FFT_nComplex; i++)
		data[i] *= norm;

	return;
}
#endif

#if defined(KGRID) && defined(FOURIERTRANSFORM)
/*
 * this is a pain because of the data
 * layout of the FFTW. For consistency we
 * store the real grid in the same way.
 */
extern void make_kGrid()
{
	double kmin;
	size_t nBytes;
	long i, j, k, i_global, iconj, jconj, kconj, idx, offset = 0;

	kmin = 2 * pi / Param.GridSize;

	nBytes = ThisTask.FFT_nComplex * sizeof(**FTGrid.kVector);

	FTGrid.kVector[0] = my_malloc(nBytes);
	FTGrid.kVector[1] = my_malloc(nBytes);
	FTGrid.kVector[2] = my_malloc(nBytes);

	if (nz % 2)
		offset = 1;

	for (i = 0; i < ThisTask.FFT_nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz_padded / 2; k++) {
				idx = (i * ny + j) * (nz_padded / 2) + k;

				i_global = i + ThisTask.FFT_x_start;

				iconj = nx - i_global;
				if (!i_global)
					iconj = 0;

				jconj = ny - j;
				if (!j)
					jconj = 0;

				kconj = nz_padded - k + offset;
				if (!k)
					kconj = 0;

				if (i_global <= nx / 2) {
					(FTGrid.kVector[0])[idx] =
					    i_global * kmin;
				} else {
					(FTGrid.kVector[0])[idx] =
					    iconj * kmin;
				}

				if (j <= ny / 2) {
					(FTGrid.kVector[1])[idx] =
					    j * kmin;
				} else {
					(FTGrid.kVector[1])[idx] =
					    jconj * kmin;
				}

				if (k <= nz_padded / 2) {
					(FTGrid.kVector[2])[idx] =
					    k * kmin;
				} else {
					(FTGrid.kVector[2])[idx] =
					    kconj * kmin;
				}
			}
		}
	}

	return;
}
#endif
