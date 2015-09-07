#include "allvars.h"

#include <hdf5.h>

#ifdef VISIT			/* Visit output */

#include <silo.h>
#define N_MESH_OPTS 7

void write_output()
{
	float *coords[3];
	long long i, nCells;
	int ndims = 3, dims[3] = { 0 };
	double *data, *data3[3], unit = 1;
	char *varnames[20];

	DBoptlist *meshOptlist = DBMakeOptlist(N_MESH_OPTS);
	DBoptlist *varOptlist = DBMakeOptlist(1);
	DBfile *dbfile = NULL;

	start_timing(CPU_OUTPUT);

	if (!ThisTask.Rank) {	/* Parallel HDF5 is a mess */

		nCells = pow(Param.NGrid, 3);

		dims[0] = Param.NGrid;
		dims[1] = Param.NGrid;
		dims[2] = Param.NGrid;

		varnames[0] = strdup("x");
		varnames[1] = strdup("y");
		varnames[2] = strdup("z");

		dbfile = DBCreate(Param.Output_File, DB_CLOBBER, DB_LOCAL,
				  "Sph2Grid", DB_PDB);

		my_assert(dbfile != NULL, "Could not create Silo file");

		coords[0] =
		    (float *) my_malloc(Param.NGrid * sizeof(*coords));
		coords[1] =
		    (float *) my_malloc(Param.NGrid * sizeof(*coords));
		coords[2] =
		    (float *) my_malloc(Param.NGrid * sizeof(*coords));

		for (i = 0; i < Param.NGrid; i++) {
			coords[0][i] = i;
			coords[1][i] = i;
			coords[2][i] = i;
		}

		DBAddOption(meshOptlist, DBOPT_DTIME,
			    (void *) &Snap.Redshift);
		DBAddOption(meshOptlist, DBOPT_XLABEL, (void *) "x");
		DBAddOption(meshOptlist, DBOPT_YLABEL, (void *) "y");
		DBAddOption(meshOptlist, DBOPT_ZLABEL, (void *) "z");
		DBAddOption(meshOptlist, DBOPT_XUNITS, (void *) "kpc");
		DBAddOption(meshOptlist, DBOPT_YUNITS, (void *) "kpc");
		DBAddOption(meshOptlist, DBOPT_ZUNITS, (void *) "kpc");

		DBPutQuadmesh(dbfile, "grid", NULL, coords, dims, ndims,
			      DB_FLOAT, DB_COLLINEAR, meshOptlist);
		free(coords[0]);
		free(coords[1]);
		free(coords[2]);

		printf("Writing MASS \n");

		unit = Unit.Mass;

		data = (double *) my_malloc(nCells * sizeof(*data));
		for (i = 0; i < nCells; i++)
			data[i] = Grid[i].Mass * unit;

		DBAddOption(varOptlist, DBOPT_UNITS, (void *) "g");
		DBPutQuadvar1(dbfile, "density", "grid", data, dims,
			      ndims, NULL, 0, DB_DOUBLE, DB_NODECENT,
			      varOptlist);

		DBClearOption(varOptlist, 0);
		my_free(data);

		printf("Writing RHO \n");

		unit = Unit.Mass / pow(Unit.Length, 3);

		data = (double *) my_malloc(nCells * sizeof(*data));
		for (i = 0; i < nCells; i++)
			data[i] = Grid[i].Rho * unit;

		DBAddOption(varOptlist, DBOPT_UNITS, (void *) "g/cm^3");
		DBPutQuadvar1(dbfile, "density", "grid", data, dims,
			      ndims, NULL, 0, DB_DOUBLE, DB_NODECENT,
			      varOptlist);

		DBClearOption(varOptlist, 0);
		my_free(data);

#ifdef VEL
		printf("VEL \n");

		unit = Unit.Vel;

		data3[0] = (double *) my_malloc(nCells * sizeof(*data3));
		data3[1] = (double *) my_malloc(nCells * sizeof(*data3));
		data3[2] = (double *) my_malloc(nCells * sizeof(*data3));
		for (i = 0; i < nCells; ++i) {
			data3[0][i] = Grid[i].Vel[0] * unit;
			data3[1][i] = Grid[i].Vel[1] * unit;
			data3[2][i] = Grid[i].Vel[2] * unit;
		}

		DBAddOption(varOptlist, DBOPT_UNITS, (void *) "cm/s");
		DBPutQuadvar(dbfile, "velocity", "grid", 3, varnames,
			     data3, dims, ndims, NULL, 0, DB_DOUBLE,
			     DB_NODECENT, varOptlist);

		DBClearOption(varOptlist, 0);
		free(data3[0]);
		free(data3[1]);
		free(data3[2]);
#endif
#ifdef MOMENTUM
		printf("Writing Momentum \n");

		unit = Unit.Vel * Unit.Mass;

		data3[0] = (double *) my_malloc(nCells * sizeof(*data3));
		data3[1] = (double *) my_malloc(nCells * sizeof(*data3));
		data3[2] = (double *) my_malloc(nCells * sizeof(*data3));
		for (i = 0; i < nCells; ++i) {
			data3[0][i] = Grid[i].Mom[0] * unit;
			data3[1][i] = Grid[i].Mom[1] * unit;
			data3[2][i] = Grid[i].Mom[2] * unit;
		}

		DBAddOption(varOptlist, DBOPT_UNITS, (void *) "g*cm/s");
		DBPutQuadvar(dbfile, "momentum", "grid", 3, varnames,
			     data3, dims, ndims, NULL, 0, DB_DOUBLE,
			     DB_NODECENT, varOptlist);

		DBClearOption(varOptlist, 0);
		free(data3[0]);
		free(data3[1]);
		free(data3[2]);
#endif
#ifdef INTENERGY
		printf("Writing U \n");

		unit = 1;

		data = (double *) my_malloc(nCells * sizeof(*data));
		for (i = 0; i < nCells; i++)
			data[i] = Grid[i].U * unit;

		DBAddOption(varOptlist, DBOPT_UNITS, (void *) "none");

		DBPutQuadvar1(dbfile, "internal_energy", "grid", data,
			      dims, ndims, NULL, 0, DB_DOUBLE, DB_NODECENT,
			      varOptlist);

		DBClearOption(varOptlist, 0);
		free(data);
#endif

#ifdef BFLD
		printf("Writing BFLD \n");

		unit = 1;

		data3[0] = (double *) my_malloc(nCells * sizeof(*data3));
		data3[1] = (double *) my_malloc(nCells * sizeof(*data3));
		data3[2] = (double *) my_malloc(nCells * sizeof(*data3));
		for (i = 0; i < nCells; ++i) {
			data3[0][i] = Grid[i].Bfld[0] * unit;
			data3[1][i] = Grid[i].Bfld[1] * unit;
			data3[2][i] = Grid[i].Bfld[2] * unit;
		}

		DBAddOption(varOptlist, DBOPT_UNITS, (void *) "G");

		DBPutQuadvar(dbfile, "magnetic_field", "grid", 3, varnames,
			     data3, dims, ndims, NULL, 0, DB_DOUBLE,
			     DB_NODECENT, varOptlist);

		DBClearOption(varOptlist, 0);
		free(data3[0]);
		free(data3[1]);
		free(data3[2]);

#endif

#ifdef VTURB
		printf("Writing VTURB \n");

		unit = Unit.Vel;

		data = (double *) my_malloc(nCells * sizeof(data));
		for (i = 0; i < nCells; i++)
			data[i] = Grid[i].VTurb * unit;

		DBAddOption(varOptlist, DBOPT_UNITS, (void *) "cm/s");

		DBPutQuadvar1(dbfile, "turbulent_velocity", "grid", data,
			      dims, ndims, NULL, 0, DB_DOUBLE, DB_NODECENT,
			      varOptlist);

		DBClearOption(varOptlist, 0);
		free(data);
#endif

#ifdef SCALARVEL
		printf("Writing SCALARVEL \n");

		unit = Unit.Vel * Unit.Vel;

		data = (double *) my_malloc(nCells * sizeof(*data));
		for (i = 0; i < nCells; i++)
			data[i] = Grid[i].VelScalar * unit;

		DBAddOption(varOptlist, DBOPT_UNITS, (void *) "cm/s");

		DBPutQuadvar1(dbfile, "Kinetic Energy v^2", "grid", data,
			      dims, ndims, NULL, 0, DB_DOUBLE, DB_NODECENT,
			      varOptlist);

		DBClearOption(varOptlist, 0);
		free(data);
#endif

		DBClose(dbfile);
	}

	stop_timing(CPU_OUTPUT);

	return;
}

#else // HDF5 Output 

void write_header(hid_t);
hid_t create_group(hid_t, char *, char *);
float *reduce_gridfield(hid_t, void *, size_t, int, int);
void write_k_grid(hid_t);
void write_h5_data(hid_t, char *, int, size_t *, int, void *, hid_t);

void write_output()
{
	hid_t fp = 0, group = 0;
	void *data;
	char group_name[MAXLINELENGTH], name[MAXLINELENGTH],
	    msg[MAXLINELENGTH];
	unsigned int flag = H5F_ACC_TRUNC;
	enum gridfields block;
	size_t dims[3] = { 0 }, nBytes;

	start_timing(CPU_OUTPUT);

	if (!ThisTask.Rank) {
		printf("Writing HDF5 file <%s> \n", Param.Output_File);

		if (Param.NoClobber)	/* handle overwrite */
			flag = H5F_ACC_EXCL;

		fp = H5Fcreate(Param.Output_File, flag, H5P_DEFAULT,
			       H5P_DEFAULT);

		sprintf(msg, "can't open file <%s>", Param.Output_File);
		my_assert(fp >= 0, msg);

		write_header(fp);
	}

	for (block = (enum gridfields) 0; block < GRID_LASTENTRY; block++) {
		set_gridfield_prop(block);

		if (!GFld.rank)
			continue;

		if (!ThisTask.Rank)
			printf("    %s \n", GFld.name);

		data = reduce_gridfield(fp, GFld.ptr, ThisTask.NDouble,
                GFld.rank, 0);

		if (!ThisTask.Rank) {
			/* Define new group */
			sprintf(group_name, "/%s", GFld.name);
			group =
			    H5Gcreate(fp, group_name, H5P_DEFAULT,
				      H5P_DEFAULT, H5P_DEFAULT);

			/* Write normal grid */
			sprintf(name, "%s/Grid", group_name);
			dims[0] = dims[1] = dims[2] = Param.NGrid;
			write_h5_data(fp, name, GFld.rank, dims, 3, data,
				      H5T_NATIVE_FLOAT);
		}
		my_free(data);

#if defined(FOURIERTRANSFORM) && !defined(NO_FFT_OUTPUT)
		/* Write FFT grids */
		data =
		    reduce_gridfield(fp, GFld.FTptr, ThisTask.FFT_nComplex,
				     GFld.rank, 1);

		if (!ThisTask.Rank && GFld.FTptr != NULL) {
			sprintf(name, "%s/FFTGrid_real", group_name);
			dims[0] = (Param.NGrid / 2 + 1);
			write_h5_data(fp, name, GFld.rank, dims, 3, data,
				      H5T_NATIVE_FLOAT);
		}
		my_free(data);

		data =
		    reduce_gridfield(fp, GFld.FTptr, ThisTask.FFT_nComplex,
				     GFld.rank, 2);

		if (!ThisTask.Rank && GFld.FTptr != NULL) {
			sprintf(name, "%s/FFTGrid_imag", group_name);
			write_h5_data(fp, name, GFld.rank, dims, 3, data,
				      H5T_NATIVE_FLOAT);
		}
		my_free(data);
#endif

#ifdef POWERSPECTRUM
		/* Write binned P(k) values */
		if (!ThisTask.Rank && GFld.Pkptr != NULL) {
			sprintf(name, "%s/Pk", group_name);

			nBytes = Param.Pk_nbins * sizeof(*GFld.Pkptr);
			data = my_malloc(nBytes);

			memcpy(data, *GFld.Pkptr, nBytes);

			dims[0] = Param.Pk_nbins;
			write_h5_data(fp, name, 1, dims, 1, data,
				      H5T_NATIVE_DOUBLE);

			my_free(data);
		}
#endif

		if (!ThisTask.Rank)
			H5Gclose(group);

		MPI_Barrier(MPI_COMM_WORLD);
	}

#if defined(KGRID) && defined(FOURIERTRANSFORM) && !defined(NO_FFT_OUTPUT)
	/* write k-grid for FFT */
	data = reduce_gridfield(fp, FTGrid.kVector, ThisTask.FFT_nComplex, 3, 3);

	if (!ThisTask.Rank) {
		printf("    KVEC \n");

		group = H5Gcreate(fp, "/KVEC", H5P_DEFAULT, H5P_DEFAULT,
			      H5P_DEFAULT);

		strcpy(name, "/KVEC/Grid");

		dims[0] = (Param.NGrid / 2. + 1);
		dims[1] = dims[2] = Param.NGrid;
		write_h5_data(fp, name, 3, dims, 3, data,
			      H5T_NATIVE_FLOAT);

		H5Gclose(group);
		my_free(data);
	}
#endif

#ifdef POWERSPECTRUM
	if (!ThisTask.Rank) { // write k-values for P(k)
		printf("    KPK \n");

		group = H5Gcreate(fp, "/KPK", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		strcpy(name, "KPK/Grid");

		dims[0] = Param.Pk_nbins;
		write_h5_data(fp, name, 1, dims, 1, Pk.k,
			      H5T_NATIVE_DOUBLE);
		H5Gclose(group);
	}
#endif

	if (!ThisTask.Rank)
		H5Fclose(fp);

	stop_timing(CPU_OUTPUT);

	MPI_Barrier(MPI_COMM_WORLD);

	return;
}

void write_header(hid_t fp)
{
	int i = 0;
	
    struct H5_header { // construct C header data 
		double Center[3];
		double GridSize;
		long NGrid;
		int Workmode;
		double Time;
		long Npart[N_part_types];
		int Cosmology;
		int Flag_Barycenter;
		char Input_File[MAXLINELENGTH];
	} head;

	for (i = 0; i < 3; i++)
		head.Center[i] = Param.Center[i];

	for (i = 0; i < N_part_types; i++)
		head.Npart[i] = Snap.Npart[i];

	strncpy(head.Input_File, Param.Input_File, MAXLINELENGTH);

	head.GridSize = Param.GridSize;
	head.NGrid = Param.NGrid;
	head.Workmode = Param.Workmode;
	head.Time = Snap.Time;
	head.Cosmology = Param.Cosmology;
	head.Flag_Barycenter = Param.Flag_Barycenter;

	/* array data types */
	hsize_t dims[1] = { N_part_types };
	hid_t lonArrayType = H5Tarray_create2(H5T_NATIVE_LONG, 1, dims);

	dims[0] = 3;
	hid_t dblArrayType = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, dims);
	hid_t charArrayType = H5Tcopy(H5T_C_S1);
	herr_t status = H5Tset_size(charArrayType, MAXLINELENGTH);

	/* header data type */
	hid_t dtype = H5Tcreate(H5T_COMPOUND, sizeof(struct H5_header));

	H5Tinsert(dtype, "Center",
		  HOFFSET(struct H5_header, Center), dblArrayType);
	H5Tinsert(dtype, "GridSize",
		  HOFFSET(struct H5_header, GridSize), H5T_NATIVE_DOUBLE);
	H5Tinsert(dtype, "NGrid",
		  HOFFSET(struct H5_header, NGrid), H5T_NATIVE_LONG);
	H5Tinsert(dtype, "Workmode",
		  HOFFSET(struct H5_header, Workmode), H5T_NATIVE_INT);
	H5Tinsert(dtype, "Time",
		  HOFFSET(struct H5_header, Time), H5T_NATIVE_DOUBLE);
	H5Tinsert(dtype, "Npart", HOFFSET(struct H5_header, Npart),
		  lonArrayType);
	H5Tinsert(dtype, "Cosmology", HOFFSET(struct H5_header, Cosmology),
		  H5T_NATIVE_INT);
	H5Tinsert(dtype, "Flag_Barycenter",
		  HOFFSET(struct H5_header, Flag_Barycenter),
		  H5T_NATIVE_INT);
	H5Tinsert(dtype, "InputFile",
		  HOFFSET(struct H5_header, Input_File), charArrayType);

	/* dataspace */
	hid_t dspace = H5Screate(H5S_SCALAR);

	/* dataset */
	hid_t dset = H5Dcreate2(fp, "HEAD", dtype, dspace,
			  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	if (!ThisTask.Rank)
		printf("    Header \n");

	status =
	    H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &head);

	/* free ids */
	H5Sclose(dspace);
	H5Dclose(dset);

	return;
}

/* data is distributed over all CPUs 
 * we combine the 3D data but store complex
 * components seperately */
float *reduce_gridfield(hid_t fp, void *data, size_t nData, int rank, int mode)
{
	size_t nDouble = 0, i, j;

	if (data == NULL)	// data doesn't exist 
		return data;
	
    float *sendbuf = NULL; // Prepare sendbuffer 

	int nSend = nData * rank;

	size_t nBytes = nSend * sizeof(*sendbuf);
	
    sendbuf = my_malloc(nBytes);

	int nRecvTotal = 0, *nRecv = NULL, *offsets = NULL; // recv book keeping

	if (!ThisTask.Rank) {  
		
        nRecv = my_malloc(ThisTask.NTask * sizeof(*nRecv));
	
        offsets = my_malloc(ThisTask.NTask * sizeof(*offsets));
	}

	MPI_Gather(&nSend, 1, MPI_INT, nRecv, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
    MPI_Reduce(&nSend, &nRecvTotal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	float *recvbuf = NULL;

	if (!ThisTask.Rank) {

		offsets[0] = 0;
		
        for (i = 1; i < ThisTask.NTask; i++)
			offsets[i] = offsets[i - 1] + nRecv[i - 1];

		recvbuf = my_malloc(nRecvTotal * sizeof(*recvbuf));
	}

	switch (mode) { // fill sendbuf
	case 0:		// real 

		nDouble = sizeof(*Grid) / sizeof(double);
		
        for (j = 0; j < rank; j++)
			for (i = 0; i < nData; i++)
				sendbuf[i * rank + j] = *(((double *) data) + i * nDouble + j);
		break;
#ifdef FOURIERTRANSFORM
	case 1:		// complex - real part 

		for (j = 0; j < rank; j++)
			for (i = 0; i < nData; i++)
				sendbuf[i * rank + j] =
				    creal((((complex double **)data)[j])[i]);
		break;
	case 2:		// complex - imag part 

		for (j = 0; j < rank; j++)
			for (i = 0; i < nData; i++)
				sendbuf[i * rank + j] = 
                    cimag((((complex double **)data)[j])[i]);
		break;
#endif
#ifdef KGRID
	case 3:		// double - grid paradigm

		for (j = 0; j < rank; j++)
			for (i = 0; i < nData; i++)
				sendbuf[i * rank + j] = (((double **) data)[j])[i];
		break;

#endif
	default:

		my_assert(0, "Mode not found");
		
        break;
	}

	MPI_Gatherv(sendbuf, nSend, MPI_FLOAT, recvbuf, nRecv, offsets, MPI_FLOAT,
            0, MPI_COMM_WORLD);

	my_free(sendbuf);
	my_free(nRecv);
	my_free(offsets);

	return recvbuf;
}

void write_h5_data(hid_t fp, char *name, int ncomp, size_t * dims, int ndims,
	      void *data, hid_t el_type)
{
	if (ThisTask.Rank)
		return;

	hid_t dtype; // first set datatype
	hsize_t elem_dims[1];

    if (ncomp == 1) {  

		dtype = el_type;
	
    } else {
	
        elem_dims[0] = ncomp;
		
        dtype = H5Tarray_create2(el_type, 1, elem_dims);
	}

	hsize_t h5dims[3] = { dims[0], dims[1], dims[2] }; // lib destroys dims[3]	 
	hid_t dspace = H5Screate_simple(ndims, h5dims, NULL); // create dataspace

	
	hid_t dset = H5Dcreate2(fp, name, dtype, dspace,  // create dataset 
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	/* write data */
	herr_t status = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

	H5Sclose(dspace); // free ids 
	H5Dclose(dset);

	return;
}

#endif
