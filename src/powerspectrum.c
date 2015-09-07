
#include "allvars.h"

#ifdef POWERSPECTRUM

#define MIN_SAMPLE_PER_BIN 50	/* approx val */

double *make_k_bins();
void bin_data_local(double *, double *, double *, long *);

void powerspectrum()
{
    int i = 0;

	const double kmin = 2 * pi / Param.GridSize;

	const int nBins = Param.Pk_nbins;

	if (!ThisTask.Rank)
		printf("Binning Powerspectrum, " "<%zu>\n", nBins);

	/* get upper bound of bins */
	double *kbins = make_k_bins();

	/* prepare buffers for local sums */
	long *local_N = my_malloc( nBins * sizeof(*local_N) );

	double *local_Pk = my_malloc( nBins * sizeof(*local_Pk) );

	/* prepare grid of k vector lengths */
    double *kgrid = my_malloc( ThisTask.FFT_nComplex * sizeof(*kgrid) );

	for (size_t i = 0; i < ThisTask.FFT_nComplex; i++)
		kgrid[i] = sqrt( p2(FTGrid.kVector[0][i]) + p2(FTGrid.kVector[1][i])
				       + p2(FTGrid.kVector[2][i]) ) / kmin;


	/* prepare MPI buffers on 0 */
	double *global_Pk = NULL;
	long  *global_N = NULL;

	if (!ThisTask.Rank) {
		size_t nBytes = nBins * sizeof(*global_Pk);
		global_Pk = my_malloc(nBytes);

		nBytes = nBins * sizeof(*global_N);
		global_N = my_malloc(nBytes);
	}

	/* bin spectrum for every field */
	for (int block = 0; block < GRID_LASTENTRY; block++) {

		set_gridfield_prop((enum gridfields) block);

		if (!GFld.rank || !GFld.Pkptr)
			continue;

		if (!ThisTask.Rank)
			printf("    %s\n", GFld.name);

		for (i = 0; i < nBins; i++)
			local_N[i] = local_Pk[i] = 0;

		bin_data_local(kbins, kgrid, local_Pk, local_N);

		if (!ThisTask.Rank)
			for (i = 0; i < nBins; i++)
				global_Pk[i] = global_N[i] = 0;

		MPI_Reduce(local_Pk, global_Pk, nBins,
			   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		MPI_Reduce(local_N, global_N, nBins,
			   MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

		if (!ThisTask.Rank) {

			size_t nBytes = nBins * sizeof(*GFld.Pkptr);
			*GFld.Pkptr = my_malloc(nBytes);

			for (i = 0; i < nBins; i++) 
				(*GFld.Pkptr)[i] = global_Pk[i] / fmax(1, global_N[i]);
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}

	my_free(kbins);
	my_free(kgrid);
	my_free(local_N);
	my_free(local_Pk);
	my_free(global_N);
	my_free(global_Pk);

	if (!ThisTask.Rank)
		printf("\n");

	return;
}

/* construct log k-bins so there is at 
 * least MIN_SAMPLE_PER_BIN points in 
 * the bin. All in units so kmin=1 */
double *make_k_bins()
{
	double *k;

	size_t nbins = Param.Pk_nbins;	/* one more than needed */

	size_t nBytes = nbins * sizeof(*k);
	k = my_malloc(nBytes);	// upper bound of bin 

	double klast = Param.NGrid/2; // in kmin
	double kfirst = 2;

	double di = log10(klast / kfirst) / (nbins - 1) ;

	for (int i = 0; i < nbins; i++) 
		k[i] = kfirst * pow(10, i * di);
	
	/* make actual bin centers */
	nBytes = Param.Pk_nbins * sizeof(*Pk.k);
	Pk.k = my_malloc(nBytes);

	double kmin = 2 * pi / Param.GridSize; 

	/* center bin position and go to physical units */
	Pk.k[0] = k[0] / 1.26 * kmin;

	for (int i = 1; i < nbins; i++) 
		Pk.k[i] = 0.5*(k[i-1] + k[i]) * kmin; 
	

	return k;
}

/* P(k) is defined as the square of the 
 * absolute of the complex value. */
void bin_data_local(double *kbins, double *kgrid, double *Pk, long *N)
{
	size_t i, j;
	double val = 0;
	
	size_t nBins = Param.Pk_nbins;

	/* Init buffers */
	for (i = 0; i < Param.Pk_nbins; i++) {
		Pk[i] = 0;
		N[i] = 0;
	}

	if (GFld.rank == 1)	/* scalars */
		for (i = 0; i < ThisTask.FFT_nComplex; i++) {
			val = p2(cabs((GFld.FTptr[0])[i]));

			for (j = 0; j < nBins; j++) {
				if (kbins[j] > kgrid[i]) {
					Pk[j] += val;
					N[j]++;
					break;
				}
			}
		}

	if (GFld.rank == 3)	/* vectors */
		for (i = 0; i < ThisTask.FFT_nComplex; i++) {
			val = p2(cabs((GFld.FTptr[0][i] + GFld.FTptr[1][i]
				       + GFld.FTptr[2][i]) / 3.0));

			for (j = 0; j < nBins; j++) {
				if (kbins[j] > kgrid[i]) {
					Pk[j] += val;
					N[j]++;
					break;
				}
			}

		}

	return;
}

#undef MIN_SAMPLE_PER_BIN

#endif				/* POWERSPECTRUM */
