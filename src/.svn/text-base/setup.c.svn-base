/* Do auxillary stuff, center particles,
 * initialise the grid
 * */
#include "allvars.h"
#include "schemes.h"

static void set_barycenter();

void center()
{
	long long ipart = 0;

	start_timing(CPU_DOMAINDECOMP);

	/* Make Center not comoving=physical */
	Param.Center[0] = (double) (Param.Center[0] * Comv2phys.Length);
	Param.Center[1] = (double) (Param.Center[1] * Comv2phys.Length);
	Param.Center[2] = (double) (Param.Center[2] * Comv2phys.Length);

	if (Param.Flag_Barycenter)
		set_barycenter();

	/* Center Snapshot on Box  */
	for (ipart = 0; ipart < ThisTask.PartTotal; ipart++) {
		P[ipart].Pos[0] =
		    (float) (P[ipart].Pos[0] - Param.Center[0]);
		P[ipart].Pos[1] =
		    (float) (P[ipart].Pos[1] - Param.Center[1]);
		P[ipart].Pos[2] =
		    (float) (P[ipart].Pos[2] - Param.Center[2]);
	}

	Param.GridSize /= (float) (Snap.Redshift + 1);	/* Make Physical */

	if (!ThisTask.Rank)
		fprintf(stdout,
			"Mesh Centered to <%6.3f,%6.3f,%6.3f> kpc physical \n"
			"Mesh Covers <%6.2f> kpc with <%6zu> cells a side \n"
			"One Cell spans <%6.2f> kpc \n\n",
			Param.Center[0], Param.Center[1], Param.Center[2],
			Param.GridSize, Param.NGrid,
			Param.GridSize / Param.NGrid);
	fflush(stdout);

	return;
}

/*Find Center of Mass
 * of the snapshot
 * */
void set_barycenter()
{
	size_t ipart;
	double send_buf[4] = { 0 }, recv_buf[4] = { 0 };

	for (ipart = 0; ipart < ThisTask.PartTotal; ipart++) {
			send_buf[0] += P[ipart].Mass * P[ipart].Pos[0];
			send_buf[1] += P[ipart].Mass * P[ipart].Pos[1];
			send_buf[2] += P[ipart].Mass * P[ipart].Pos[2];
            send_buf[3] += P[ipart].Mass;
	}

	MPI_Allreduce(send_buf, recv_buf, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double mtot = recv_buf[3];
    my_assert(mtot, "Failed to compute total mass");

	recv_buf[0] /= mtot;
	recv_buf[1] /= mtot;
	recv_buf[2] /= mtot;

	if (Param.Flag_Barycenter == 2) {
		if (!ThisTask.Rank)
			printf("Setting Center relative to Barycenter \n");

		Param.Center[0] += recv_buf[0];
		Param.Center[1] += recv_buf[1];
		Param.Center[2] += recv_buf[2];
	} else {
		if (!ThisTask.Rank)
			printf("Setting center to Barycenter \n");

		Param.Center[0] = recv_buf[0];
		Param.Center[1] = recv_buf[1];
		Param.Center[2] = recv_buf[2];
	}

	return;
}

/* fill GFld from block gfield 
 * we handle grid access and reduction with this */
void set_gridfield_prop(enum gridfields field)
{
	char msg[MAXLINELENGTH];

	GFld.rank = 0;		/* == not present */
	strncpy(GFld.name, "\0", MAXLINELENGTH);

	switch (field) {
	case GRID_MASS:
#ifdef MASS
		GFld.ptr = &Grid[0].Mass;
		GFld.FTptr = &FTGrid.Mass;
		GFld.Pkptr = &Pk.Mass;
		GFld.rank = 1;
		strncpy(GFld.name, "MASS", 4);
#endif
		break;
	case GRID_RHO:
#ifdef RHO
		GFld.ptr = &Grid[0].Rho;
		GFld.FTptr = &FTGrid.Rho;
		GFld.Pkptr = &Pk.Rho;
		GFld.rank = 1;
		strncpy(GFld.name, "RHO", 4);
#endif
		break;
	case GRID_VEL:
#ifdef VEL
		GFld.ptr = &Grid[0].Vel[0];
		GFld.FTptr = &FTGrid.Vel[0];
		GFld.Pkptr = &Pk.Vel;
		GFld.rank = 3;
		strncpy(GFld.name, "VEL", 4);
#endif
		break;
	case GRID_MOM:
#ifdef MOMENTUM
		GFld.ptr = &Grid[0].Mom[0];
		GFld.FTptr = &FTGrid.Mom[0];
		GFld.Pkptr = &Pk.Mom;
		GFld.rank = 3;
		strncpy(GFld.name, "MOM", 4);
#endif
		break;
	case GRID_U:
#ifdef INTENERGY
		GFld.ptr = &Grid[0].U;
		GFld.FTptr = &FTGrid.U;
		GFld.Pkptr = &Pk.U;
		GFld.rank = 1;
		strncpy(GFld.name, "U", 4);
#endif
		break;
	case GRID_BFLD:
#ifdef BFLD
		GFld.ptr = &Grid[0].Bfld[0];
		GFld.FTptr = &FTGrid.Bfld[0];
		GFld.Pkptr = &Pk.Bfld;
		GFld.rank = 3;
		strncpy(GFld.name, "BFLD", 4);
#endif
		break;
	case GRID_B:
#ifdef SCALAR_BFLD
		GFld.ptr = &Grid[0].B;
		GFld.FTptr = &FTGrid.B;
		GFld.Pkptr = &Pk.B;
		GFld.rank = 1;
		strncpy(GFld.name, "B", 4);
#endif
		break;
    case GRID_VTURB:
#ifdef VTURB
		GFld.ptr = &Grid[0].VTurb[0];
		GFld.FTptr = &FTGrid.VTurb;
		GFld.Pkptr = &Pk.VTurb;
		GFld.rank = 3;
		strncpy(GFld.name, "VTURB", 4);
#endif
		break;
	case GRID_NPART:
		GFld.ptr = &Grid[0].Npart;
		GFld.FTptr = NULL;
		GFld.Pkptr = NULL;
		GFld.rank = 1;
		strncpy(GFld.name, "NPRT", 4);
		break;
	case GRID_DENSVEL:
#ifdef DENSVEL
		GFld.ptr = &Grid[0].DensVel[0];
		GFld.FTptr = &FTGrid.DensVel[0];
		GFld.Pkptr = &Pk.DensVel;
		GFld.rank = 3;
		strncpy(GFld.name, "DWVEL", 4);
#endif
		break;
	case GRID_SCALARVEL:
#ifdef SCALARVEL
		GFld.ptr = &Grid[0].VelScalar;
		GFld.FTptr = &FTGrid.VelScalar;
		GFld.Pkptr = &Pk.VelScalar;
		GFld.rank = 1;
		strncpy(GFld.name, "SCAV", 4);
#endif
		break;
	case GRID_VELDISPERSION:
#ifdef VELDISPERSION
		GFld.ptr = &Grid[0].Sigma[0];
		GFld.FTptr = NULL;
		GFld.Pkptr = NULL;
		GFld.rank = 3;
		strncpy(GFld.name, "VSIG", 4);
#endif
		break;
	case GRID_VELFILTERED:
#ifdef VELFILTERED
		GFld.ptr = &Grid[0].VelFiltered[0];
		GFld.FTptr = &FTGrid.VelFiltered[0];
		GFld.Pkptr = &Pk.VelFiltered;
		GFld.rank = 3;
		strncpy(GFld.name, "VFIL", 4);
#endif
		break;
	default:
		sprintf(msg, "Gridfield %d not found", (int) field);
		my_assert(0, msg);
		break;
	}

	return;
}

/* set all grid fields to 0 */
void init_grid()
{
    const size_t nBytes = NGrid3 * sizeof(*Grid);
	Grid = my_malloc(nBytes);

    if (!ThisTask.Rank)
        printf("Initialising Grid (%g MB) \n", round(nBytes/1024./1024.0));

    memset(Grid, 0, nBytes);

	MPI_Barrier(MPI_COMM_WORLD);
	
    return;
}
