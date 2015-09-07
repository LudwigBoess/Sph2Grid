/* This is SPH2Grid - bin particle to a grid */

#include "allvars.h"
#include "input.h"
#include "tree.h"
#include "ngb.h"
#include "schemes.h"

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask.Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &ThisTask.NTask);

#ifdef OPENMP
#pragma omp parallel
	{
		Omp.ThreadID = omp_get_thread_num();
		Omp.NThreads = omp_get_num_threads();
	}
#endif

	if (!ThisTask.Rank) {
		printf("*** This is SPH2Grid Version 1.2 ***\n\n");

		printf("Running on <%d> Nodes \n", ThisTask.NTask);
#ifdef OPENMP
#pragma omp parallel
		if (!Omp.ThreadID)
			printf("Running on <%d> Threads per Node\n\n",
			       Omp.NThreads);
#endif

		print_compile_time_settings();

		my_assert(argc == 2, "One parameter file please");

	}

	init_timing();

	read_param_file(argv[1]);

	set_cosmology(Param.Cosmology);

	set_units();

	read_snapshot(Param.Input_File);

	center();

	domain_decomposition();

#ifdef TREE			/* serial tree ... */
	init_tree();

	build_tree();
#endif

#if defined(HSMLFIND) && defined(TREE)
	find_hsml();
#endif

	if (!ThisTask.Rank)
		printf("Selecting Distribution Scheme <%s>\n%s\n",
		       WORKMODE, WMODE_MSG);
	WMODE_INIT;

	init_grid();

	fill_grid();

	distribute_grid();

	test_grid();

#ifdef FOURIERTRANSFORM
	fft_grid();
#endif

#if defined(FOURIERTRANSFORM) && defined(KGRID)
	make_kGrid();
#endif

#if defined(FOURIERTRANSFORM) && defined(KGRID) && defined(POWERSPECTRUM)
	powerspectrum();
#endif

#ifdef TURBVEL
	local_turbulence();
#endif

	write_output();

	finish_timing();

	MPI_Finalize();

	return (EXIT_SUCCESS);	/* Finish him of - Fatality */
}
