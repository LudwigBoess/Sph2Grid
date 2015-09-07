/* Remove uneeded Particles from Memory.
 * Distribute the rest evenly over all CPUs
 * */

#include "allvars.h"
#include "input.h"
#include <gsl/gsl_sort.h>

static bool particle_in_mesh(long long);
static double load_above_mean();
static void xChange_particles(int, double, int, double);
static double metric(long long);
static void sort_particles();

void domain_decomposition()
{
	int i = 0, task = 0, cnt = 0, xChangeTask = 0;

	long long src = 0, target = 0;
	long long nPart[N_part_types] = { 0 };

	double total_load = 0, imbalance = 0, total_imbalance = 0;

	MPI_Barrier(MPI_COMM_WORLD);

	if (!ThisTask.Rank)
		printf("Starting Domain Decomposition\n");

	/* Remove abundant particles from memory */
	for (src = 0; src < ThisTask.PartTotal; src++) {
		if (particle_in_mesh(src)) {
			memmove(&P[target].Pos[0], &P[src].Pos[0],
				sizeof(*P));

			if (P[src].Type == 0)	// Treat Gas particles
				memmove(&Gas[target], &Gas[src],
					sizeof(*Gas));

			target++;
			nPart[P[src].Type]++;
		}
	}

	for (i = 0; i < N_part_types; i++) 
		nPart[i] = nPart[i] - ThisTask.Npart[i];

    Reallocate_P(target - ThisTask.PartTotal, nPart, -1);

	MPI_Allreduce(&ThisTask.PartTotal, &Snap.PartTotal, 1,
		      MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(ThisTask.Npart, Snap.Npart, N_part_types,
		      MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

	if (!ThisTask.Rank)
		printf("\nTotal Number of Particles on Mesh : %lli \n"
		       "Number of Gas   Particles on Mesh : %lli \n"
		       "Number of DM    Particles on Mesh : %lli \n"
		       "Number of Disk  Particles on Mesh : %lli \n"
		       "Number of Bulge Particles on Mesh : %lli \n"
		       "Number of Star  Particles on Mesh : %lli \n"
		       "Number of Bndry Particles on Mesh : %lli \n\n",
		       Snap.PartTotal, Snap.Npart[0], Snap.Npart[1],
		       Snap.Npart[2], Snap.Npart[3], Snap.Npart[4],
		       Snap.Npart[5]);
	fflush(stdout);

	my_assert(Snap.PartTotal > 0, "No Particles on Mesh");

    sort_particles(); // put the large ones first
	
	/* Distribute particles evenly */
    double *excess = my_malloc(ThisTask.NTask * sizeof(*excess));

	size_t *sort = my_malloc(ThisTask.NTask * sizeof(*sort));

	for (;;) {

		double myExcess = load_above_mean(&total_load);	// returns total_load

		MPI_Allgather(&myExcess, 1, MPI_DOUBLE,
			      excess, 1, MPI_DOUBLE, MPI_COMM_WORLD);

		for (imbalance = task = 0; task < ThisTask.NTask; task++)
			imbalance += fabs(excess[task]) / total_load;

		if (imbalance < 0.1 || cnt == 50)
			break;

		gsl_sort_index(sort, excess, 1, ThisTask.NTask);

		for (i = 0; i < ThisTask.NTask; i++) 
			if (sort[i] == ThisTask.Rank)
				xChangeTask = sort[ThisTask.NTask - i - 1];

		if (excess[ThisTask.Rank] * excess[xChangeTask] < 0) 
			xChange_particles(xChangeTask, excess[xChangeTask],
					  ThisTask.Rank, excess[ThisTask.Rank]);

		cnt++;
	}

	my_free(excess);
	my_free(sort);

	MPI_Reduce(&imbalance, &total_imbalance, 1,
		   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	if (!ThisTask.Rank) 
		printf("\nAverage imbalance of <%g> after <%d> iterations \n\n",
		     total_imbalance/ThisTask.NTask, cnt);

	stop_timing(CPU_DOMAINDECOMP);

	return;
}

/* check if particle overlaps with grid */
static bool particle_in_mesh(long long ipart)
{
	const double halfGridSize = 0.5 * Param.GridSize;
	const float boxsize = (float) Snap.Boxsize;

	if (P[ipart].Type != 0)	/* Gas only */
		return false;

#ifdef PERIODIC			// Periodic images have to be considered
	for (int k = 0; k < 8; k++) {	/* Mirror Particle */
		float px = !(k & 0x1) ?
		    	P[ipart].Pos[0] : (P[ipart].Pos[0] > 0 ?
				P[ipart].Pos[0] - boxsize : 
				P[ipart].Pos[0] + boxsize);
		float py = !(k & 0x2) ? 
			P[ipart].Pos[1] : (P[ipart].Pos[1] > 0 ? 
				P[ipart].Pos[1] - boxsize :
				P[ipart].Pos[1] + boxsize);
		float pz = !(k & 0x4) ? 
			P[ipart].Pos[2] : (P[ipart].Pos[2] > 0 ? 
				P[ipart].Pos[2] - boxsize :
				P[ipart].Pos[2] + boxsize);
#else
	float px = P[ipart].Pos[0];
	float py = P[ipart].Pos[1];
	float pz = P[ipart].Pos[2];
#endif // PERIODIC

	if (px > halfGridSize || px < -halfGridSize
	 || py > halfGridSize || py < -halfGridSize
	 || pz > halfGridSize || pz < -halfGridSize) 
#ifdef PERIODIC
		continue;
#else
        return false;
#endif // PERIODIC
    else
        return true;

#ifdef PERIODIC
    } // for (int k = 0; k < 8; k++)
#endif

    return false;
}

/* 
 * Determine excess computational 
 * load present on this CPU
 */
static double load_above_mean(double *return_ptr)
{
	long long ipart, task;
	double *load, myLoad = 0, total = 0, mean = 0, excess = 0;

	for (ipart = 0; ipart < ThisTask.Npart[0]; ipart++)
		myLoad += metric(ipart);

	load = my_malloc(ThisTask.NTask * sizeof(*load));
	MPI_Allgather(&myLoad, 1, MPI_DOUBLE, load, 1, MPI_DOUBLE,
		      MPI_COMM_WORLD);

	for (task = 0; task < ThisTask.NTask; task++)
		total += load[task];

	mean = total / ThisTask.NTask;

	excess = load[ThisTask.Rank] - mean;

    free(load);

	*return_ptr = total;	/* returns total load and excess load */

	return excess;
}

static void xChange_particles(int yourRank, double yourExcess, int myRank,
	double myExcess)
{
	long long partTotal = 0, nPart[N_part_types] = { 0 };
	int tag = 0;
	MPI_Status status;

	if (myExcess > 0) {	// Send particles 

		ptrdiff_t ipart = ThisTask.PartTotal - 1;
		
	    double sendLoad = 0;

        while (sendLoad < myExcess && sendLoad < fabs(yourExcess)) {
		
            sendLoad += metric(ipart);
			
            partTotal++;
			
            nPart[P[ipart].Type]++;
			
            ipart--;
		}

		tag = myRank;

		MPI_Ssend(&partTotal, 1, MPI_LONG_LONG, yourRank, tag,
			  MPI_COMM_WORLD);
		
        MPI_Ssend(nPart, N_part_types, MPI_LONG_LONG, yourRank,
			  tag, MPI_COMM_WORLD);

		MPI_Ssend(&(P[ThisTask.PartTotal - partTotal]),
			  partTotal * sizeof(*P), MPI_BYTE,
			  yourRank, tag, MPI_COMM_WORLD);
		
        MPI_Ssend(&(Gas[ThisTask.Npart[0] - nPart[0]]),
			  nPart[0] * sizeof(*Gas), MPI_BYTE,
			  yourRank, tag, MPI_COMM_WORLD);

		Reallocate_P(partTotal, nPart, -1);
	
    } else { // Receive particles 

		tag = yourRank;

		MPI_Recv(&partTotal, 1, MPI_LONG_LONG, yourRank,
			 tag, MPI_COMM_WORLD, &status);
		
        MPI_Recv(nPart, N_part_types, MPI_LONG_LONG,
			 yourRank, tag, MPI_COMM_WORLD, &status);

		Reallocate_P(partTotal, nPart, +1);

		MPI_Recv(&(P[ThisTask.PartTotal - partTotal]),
			 partTotal * sizeof(*P), MPI_BYTE,
			 yourRank, tag, MPI_COMM_WORLD, &status);
		
        MPI_Recv(&(Gas[ThisTask.Npart[0] - nPart[0]]),
			 nPart[0] * sizeof(*Gas), MPI_BYTE,
			 yourRank, tag, MPI_COMM_WORLD, &status);
	}

	return;
}

static void sort_particles() // acts on the last nRandom particles
{
    size_t ipart;

    const size_t nPart = ThisTask.Npart[0];

    double *cost = my_malloc(nPart * sizeof(double));
    size_t *idx = my_malloc(nPart * sizeof(size_t));

    for (ipart=0; ipart<nPart; ipart++) 
        cost[ipart] = 1.0 / metric(ipart);
  
	gsl_sort_index(idx, cost, 1, nPart);

    for (ipart=0; ipart<nPart; ipart++) {
        if (idx[ipart] == ipart)
            continue;

        struct Particle_Data Psrc = P[ipart];
        struct Gas_Data GasSrc = Gas[ipart];
        size_t dest = idx[ipart];

        for (;;) {
            struct Particle_Data Pnext = P[dest];
            struct Gas_Data GasNext = Gas[dest];
            size_t idxNext = idx[dest];

            memcpy(&P[dest], &Psrc, sizeof(*P));
            memcpy(&Gas[dest], &GasSrc, sizeof(*Gas));
            idx[dest] = dest;

            if (dest == ipart)
                break;

            memcpy(&Psrc, &Pnext, sizeof(*P));
            memcpy(&GasSrc, &GasNext, sizeof(*Gas));
            dest = idxNext;

        }
    }

    return;
}

/* 
 * Return a measure of the computational
 * load carried by particle ipart
 */

static double metric(long long ipart)
{
    float rho_avg = Snap.Npart[0]/Snap.Boxsize;
#ifdef SPH
	return Gas[ipart].Rho / rho_avg;
#else
    return 1; // constant binning cost
#endif
}
