#include "allvars.h"
#include "input.h"
#include <limits.h>

#define int4bytes int

#define SKIPF77HEAD  {if(my_fread(&blksize,sizeof(int),1,fd)){swap_Nbyte((char*)&blksize,1,4);}}

#define UINT sizeof(unsigned int)
#define INT sizeof(long)
#define FLOAT sizeof(float)
#define DBL sizeof(double)

int4bytes blksize, swap = 0;

static long n_files;

/* 
 * Main Reading :
 * This routine tries first to put one snapshot per task and read
 * in N_IOTasks groups to save communication.
 * Then N_IOTasks Masters distribute to groups
 * All files left are then read at the same time.
 */
extern void read_snapshot(char *filename)
{
	long rest_files, groupMaster, groupLast, groupSize, group,
	    file_number;
	char buf[MAXLINELENGTH];

	start_timing(CPU_READIN);

	Snap.Have_Arepo = 0;

	n_files = rest_files = find_files(filename);

	if (Param.N_IOTasks > n_files) {
		Param.N_IOTasks = n_files;
		if (!ThisTask.Rank)
			printf("\nReducing # Read Tasks to # Files \n\n");
	}

	groupSize = ThisTask.NTask / Param.N_IOTasks;
	if (ThisTask.NTask % Param.N_IOTasks)
		groupSize++;

	groupMaster = (ThisTask.Rank / groupSize) * groupSize;

	groupLast = groupMaster + groupSize - 1;
	if (groupLast > ThisTask.NTask - 1)
		groupLast = ThisTask.NTask - 1;

	while (rest_files > 0) {
		if (n_files == 1) {	/* no parallel read */
			sprintf(buf, "%s", filename);

			read_file(buf, 0, ThisTask.NTask - 1);

			MPI_Barrier(MPI_COMM_WORLD);

			rest_files -= 1;
		} else if (rest_files >= ThisTask.NTask) {	/* parallel no communic. */
			file_number =
			    ThisTask.Rank + (rest_files - ThisTask.NTask);

			sprintf(buf, "%s.%li", filename, file_number);

			for (group = 0; group < groupSize; group++) {

				if (ThisTask.Rank == (groupMaster + group))
					read_file(buf, ThisTask.Rank,
						  ThisTask.Rank);

				MPI_Barrier(MPI_COMM_WORLD);
			}

			rest_files -= ThisTask.NTask;
		} else if (rest_files >= Param.N_IOTasks) {	/* parallel & distrib. */
			file_number =
			    groupMaster / groupSize + (rest_files -
						       Param.N_IOTasks);

			sprintf(buf, "%s.%li", filename, file_number);

			read_file(buf, groupMaster, groupLast);

			MPI_Barrier(MPI_COMM_WORLD);

			rest_files -= Param.N_IOTasks;
		} else {	/* parallel & distrib. of all files left  */
			groupSize = ThisTask.NTask / rest_files;

			if (ThisTask.NTask % Param.N_IOTasks)
				groupSize++;

			groupMaster =
			    (ThisTask.Rank / groupSize) * groupSize;

			file_number = groupMaster / groupSize;

			sprintf(buf, "%s.%li", filename, file_number);

			read_file(buf, groupMaster, groupLast);

			MPI_Barrier(MPI_COMM_WORLD);

			rest_files -= groupSize;
		}
	}

	if (Snap.Have_Arepo && !ThisTask.Rank)
		fprintf(stdout, "Read an AREPO Snapshot ! \n\n");

	stop_timing(CPU_READIN);

	return;
}

/* Reads and Distributes a file 
 * */
void read_file(char *filename, int ReadTask, int LastTask)
{
	long i = 0, j = 0, task = 0;
	long target, src, nTask, blockExist;
	long long nRecv[N_part_types] = { 0 }, nReadTot = 0;
	long long nSend[N_part_types] = { 0 }, ntot = 0;
	char *comm_buf = NULL;
	size_t nBytes, byteOffset, partOffset, bufOffset;
	FILE *fd = NULL;
	enum iofields blocknr;
	int tag = ReadTask;
	MPI_Status status;

	if (ThisTask.Rank == ReadTask) {
		printf("Reading file <%s> on Task <%i-%i> \n",
                filename, ReadTask, LastTask); 

		fd = fopen(filename, "r");
		ntot = read_gadget_head(fd);

		for (task = ReadTask + 1; task <= LastTask; task++)
			MPI_Ssend(&Header, sizeof(Header), MPI_BYTE, task,
				  tag, MPI_COMM_WORLD);

	} else 
		MPI_Recv(&Header, sizeof(Header), MPI_BYTE, ReadTask, tag,
			 MPI_COMM_WORLD, &status);
    

	if (!ThisTask.Rank)
		printf("Total Number of Particles in file : %lli \n"
		       "Total Number of Gas   Particles in file : %lli \n"
		       "Total Number of DM    Particles in file : %lli \n"
		       "Total Number of Disk  Particles in file : %lli \n"
		       "Total Number of Bulge Particles in file : %lli \n"
		       "Total Number of Stars Particles in file : %lli \n"
		       "Total Number of Bndry Particles in file : %lli \n\n",
		       ntot, Header.Npart[0], Header.Npart[1], Header.Npart[2],
		       Header.Npart[3], Header.Npart[4], Header.Npart[5]);

	if (Snap.Boxsize == 0 && !ThisTask.Rank)	/* Only once */
			printf("Snapshot Properties \n"
			       "   Boxsize : %g \n"
			       "   Redshift : %g\n"
			       "   Time : %g\n\n",
			       Header.BoxSize, Header.Redshift, Header.Time);

	/* Set Snapshot Properties */
	Snap.Boxsize = Header.BoxSize;
	Snap.Redshift = Header.Redshift;
	Snap.Time = Header.Time;

	for (i = 0; i < N_part_types; i++) {
		Snap.Npart[i] = Header.Nall[i];
		Snap.Masstab[i] = Header.Mass[i] / Cosmo.h;
	}

	/* Set Comoving units for Gadget */
	Comv2phys.Length = 1 / (1 + Snap.Redshift) / Cosmo.h;
	Comv2phys.Mass = 1 / Cosmo.h;
	Comv2phys.Vel = sqrt(1. / (1 + Snap.Redshift));

	/* Determine particle distribution over CPUs */
	nTask = LastTask - ReadTask + 1;
	for (i = 0; i < N_part_types; i++) {
		for (j = ThisTask.Rank - ReadTask; j < Header.Npart[i];
		     j += nTask) {
			nRecv[i]++;
			nReadTot++;
		}
	}

	/* make space for new particles */
	Reallocate_P(nReadTot, nRecv, +1);

	/* Shift collisionless particles if multiple files */
	src = ThisTask.Npart[0] - nRecv[0];
	target = src + nReadTot;
	nBytes = (ThisTask.PartTotal - nReadTot - src) * sizeof(*P);
	memmove(&P[target].Pos[0], &P[src].Pos[0], nBytes);

	/* Read blocks  */
	for (blocknr = 0; blocknr < IO_LASTENTRY; blocknr++) {
		set_block_prop(blocknr);

		byteOffset = 0;

		if (ThisTask.Rank == ReadTask) {

			blockExist = find_block(fd, Block.Label);

			if (blockExist) {
				comm_buf =
				    my_malloc(Block.Ntot *
					      Block.Bytes_per_element);

				blockExist =
				    read_gadget_block(comm_buf,
						      Block.Label, fd,
						      Block.Data_type);
			}

			for (task = ReadTask + 1; task < LastTask + 1;
			     task++)
				MPI_Ssend(&blockExist, 1, MPI_LONG, task,
					  tag, MPI_COMM_WORLD);
		} else {
			MPI_Recv(&blockExist, 1, MPI_LONG, ReadTask, tag,
				 MPI_COMM_WORLD, &status);
        }

		/* Distribute data */
		for (i = 0; i < N_part_types; i++) {

			if (!blockExist)
				continue;	/* nothing to do */

			if (!Block.Npart[i])
				continue;	/* block is not in these types */

			if (ThisTask.Rank == ReadTask) {

				byteOffset +=
				    nRecv[i] * Block.Bytes_per_element;

				for (task = ReadTask + 1;
				     task < LastTask + 1; task++) {
					MPI_Recv(nSend, N_part_types,
						 MPI_LONG_LONG, task, tag,
						 MPI_COMM_WORLD, &status);

					nBytes =
					    nSend[i] *
					    Block.Bytes_per_element;

					MPI_Ssend(comm_buf + byteOffset,
						  nBytes, MPI_BYTE, task,
						  tag, MPI_COMM_WORLD);

					byteOffset += nBytes;
				}
			} else {	/* ThisTask.Rank == RecvTask */

				if (!i){
                    nBytes = Block.Bytes_per_element*nReadTot;
					comm_buf = my_malloc(nBytes);
                }

				MPI_Ssend(nRecv, N_part_types,
					  MPI_LONG_LONG, ReadTask, tag,
					  MPI_COMM_WORLD);

				nBytes = nRecv[i] * Block.Bytes_per_element;

				MPI_Recv(comm_buf + byteOffset, nBytes,
					 MPI_BYTE, ReadTask, tag,
					 MPI_COMM_WORLD, &status);

				byteOffset += nBytes;
			}
		}

		/* transfer data to particles */
		partOffset = ThisTask.Npart[0] - nRecv[0];
		bufOffset = 0;
		for (i = 0; i < N_part_types; i++) {
			if (Block.Npart[i] && blockExist) {

				for (j = 0; j < nRecv[i]; j++)
					empty_comm_buffer(blocknr,
							  comm_buf, partOffset + j,
							  (bufOffset + j) *
							  Block.Val_per_element);

				if (blocknr == IO_ID)	/* set type */
					for (j = 0; j < nRecv[i]; j++)
						P[partOffset + j].Type = i;

			/* no MASS block so masses from Header */
			} else if (blocknr == IO_MASS) {

				for (j = 0; j < nRecv[i]; j++)
					P[partOffset + j].Mass =
					    (float) Snap.Masstab[i];
            }


			if (Block.Npart[i]) {
				if (ThisTask.Rank == ReadTask)
					bufOffset += Block.Npart[i];
				else
					bufOffset += nRecv[i];

				partOffset += nRecv[i];
			}
		}
		if (blockExist)
			free(comm_buf);
	}

	/* treat AREPO VOL -> HSML */
	if (Snap.Have_Arepo)
		for (j = 0; j < nRecv[0]; j++)
			Gas[j].Hsml =
			    (float) 4 *pow(Gas[j].Hsml / fourpithirds,
					   1.0 / 3.0);

	if (ThisTask.Rank == ReadTask)
		fclose(fd);

	return;
}

/*Determine number of files to read
 * */
int find_files(char *filename)
{
	char buf[MAXLINELENGTH], msg[MAXLINELENGTH];
	int n_files = 0;
	FILE *fd = NULL;

	if (!(fd = fopen(filename, "r"))) {
		for (;;) {
			sprintf(buf, "%s.%i", filename, n_files);
			if (!(fd = fopen(buf, "r"))) 
				break;
			
			fclose(fd);
			n_files++;

			my_assert(n_files < 1000, "Too many files");
		}

		sprintf(msg, "Can't open file <%s> or <%s>", filename, buf);
		my_assert(n_files, msg);

		if (!ThisTask.Rank)
			printf("\n Found <%i> file(s) ! \n\n", n_files);

	} else 
		n_files = 1;

	return n_files;
}

/* Basic routine to read data from a file 
 * */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
	size_t nRead = fread(ptr, size, nmemb, stream);

	if (nRead != nmemb) {
		my_assert(feof(stream), "I/O error (fread) !");

		nRead = 0;	/*EOF reached */
	}

	return nRead;
}

/* Routine to swap ENDIAN
 * */
void swap_Nbyte(char *data, int n, int m)
{
	int i, j;
	char old_data[16];

	if (swap > 0) {
		for (j = 0; j < n; j++) {
			memcpy(&old_data[0], &data[j * m], m);
			for (i = 0; i < m; i++) {
				data[j * m + i] = old_data[m - i - 1];
			}
		}
	}

	return;
}

int find_block(FILE * fd, char *label)
{
	int4bytes blocksize = 0;
	char blocklabel[5] = { "    " };

	rewind(fd);

	while (!feof(fd) && blocksize == 0) {
		SKIPF77HEAD;
		if (blksize == 134217728) {
			printf("Enabling ENDIAN swapping !\n");
			swap = 1 - swap;
			swap_Nbyte((char *) &blksize, 1, 4);
		}

		my_assert(blksize == 8, "Incorrect file format");

		if (my_fread(blocklabel, 4 * sizeof(char), 1, fd)) {
			my_fread(&blocksize, sizeof(int4bytes), 1, fd);
			swap_Nbyte((char *) &blocksize, 1, 4);
			SKIPF77HEAD;
			if (strcmp(label, blocklabel) != 0) {
				fseek(fd, blocksize, 1);
				blocksize = 0;
			}
		} else {
			blocksize = 8;
			break;
		}
	}
	return (blocksize - 8);
}

/* Read the header information, returns total no of particles in file
 * */
int read_gadget_head(FILE * fd)
{
	int blocksize, dummysize, i;
	unsigned int Npart[N_part_types], Nall[N_part_types],
	    NallHW[N_part_types];
	long long ntot = 0;

	blocksize = find_block(fd, "HEAD");

	my_assert(blocksize > 0, "Header not found");

	dummysize = blocksize - 2 * N_part_types * sizeof(int) -
	    4 * sizeof(long) - 12 * sizeof(double);
	SKIPF77HEAD;

	my_fread(Npart, N_part_types * sizeof(unsigned int), 1, fd);
	swap_Nbyte((char *) Npart, N_part_types, 4);

	my_fread(Header.Mass, N_part_types * sizeof(double), 1, fd);
	swap_Nbyte((char *) Header.Mass, N_part_types, 8);

	my_fread((void *) &Header.Time, sizeof(double), 1, fd);
	swap_Nbyte((char *) &Header.Time, 1, 8);

	my_fread((void *) &Header.Redshift, sizeof(double), 1, fd);
	swap_Nbyte((char *) &Header.Redshift, 1, 8);

	my_fread((void *) &Header.FlagSfr, sizeof(int), 1, fd);
	swap_Nbyte((char *) &Header.FlagSfr, 1, 4);

	my_fread((void *) &Header.FlagFeedback, sizeof(int), 1, fd);
	swap_Nbyte((char *) &Header.FlagFeedback, 1, 4);

	my_fread(Nall, N_part_types * sizeof(unsigned int), 1, fd);
	swap_Nbyte((char *) Nall, N_part_types, 4);

	my_fread((void *) &Header.FlagCooling, sizeof(int), 1, fd);
	swap_Nbyte((char *) &Header.FlagCooling, 1, 4);

	my_fread((void *) &Header.NumFiles, sizeof(int), 1, fd);
	swap_Nbyte((char *) &Header.NumFiles, 1, 4);

	my_fread((void *) &Header.BoxSize, sizeof(double), 1, fd);
	swap_Nbyte((char *) &Header.BoxSize, 1, 8);

	my_fread((void *) &Header.Omega0, sizeof(double), 1, fd);
	swap_Nbyte((char *) &Header.Omega0, 1, 8);

	my_fread((void *) &Header.OmegaLambda, sizeof(double), 1, fd);
	swap_Nbyte((char *) &Header.OmegaLambda, 1, 8);

	my_fread((void *) &Header.HubbleParam, sizeof(double), 1, fd);
	swap_Nbyte((char *) &Header.HubbleParam, 1, 8);

	my_fread((void *) &Header.FlagAge, sizeof(int), 1, fd);
	swap_Nbyte((char *) &Header.FlagAge, 1, 8);

	my_fread((void *) &Header.FlagMetals, sizeof(int), 1, fd);
	swap_Nbyte((char *) &Header.FlagMetals, 1, 8);

	my_fread((void *) NallHW, N_part_types * sizeof(unsigned int), 1,
		 fd);
	swap_Nbyte((char *) NallHW, N_part_types, 4);

	if (NallHW[0] != 0)
		printf("Nall HW not well tested !! ");

	fseek(fd, dummysize, 1);
	SKIPF77HEAD;

	for (i = 0; i < N_part_types; i++) {	/* HighWord */
		Header.Npart[i] =
		    (long long) (Npart[i] +
				 (((long long) NallHW[i]) << 32));
		Header.Nall[i] =
		    (long long) (Nall[i] +
				 (((long long) NallHW[i]) << 32));
		ntot += Header.Npart[i];
	}

	return ntot;
}

int
read_gadget_block(void *data, char *label, FILE * fd, size_t sizeof_type)
{
	int blocksize = 0;

	blocksize = find_block(fd, label);
	if (blocksize <= 0) {
		return (blocksize);
	} else {
		SKIPF77HEAD;

		my_fread(data, blocksize, 1, fd);
		
        swap_Nbyte((char *) data, blocksize / sizeof_type, 4);
		
        SKIPF77HEAD;
	}

	return blocksize;
}

/* Set Block characteristics 
 * "You are all different" - 
 * "We are all different" - 
 * "I'm not !" 
 * 				(life of brian) */
void set_block_prop(enum iofields blocknr)
{
	int i;

	for (i = 0; i < N_part_types; i++)
		Block.Npart[i] = 0;

	switch (blocknr) {
	case IO_POS:
		Block.Label = "POS ";	/* Has to be 4 Letters */
		Block.Name = "Coordinates";
		for (i = 0; i < N_part_types; i++)
			Block.Npart[i] = Header.Npart[i];
		Block.Val_per_element = 3;	/* 3=vector, 1=scalar */
		Block.Data_type = FLOAT;	/* Length in bytes */
		Block.Rmv_comoving = Comv2phys.Length;
		break;
	case IO_VEL:
		Block.Label = "VEL ";
		Block.Name = "Velocities";
		for (i = 0; i < N_part_types; i++)
			Block.Npart[i] = Header.Npart[i];
		Block.Val_per_element = 3;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Vel;
		break;
	case IO_ID:
		Block.Label = "ID  ";
		Block.Name = "Particle IDs";
		for (i = 0; i < N_part_types; i++)
			Block.Npart[i] = Header.Npart[i];
		Block.Val_per_element = 1;
		Block.Data_type = UINT;
		Block.Rmv_comoving = 1;
		break;
	case IO_U:
		Block.Label = "U   ";
		Block.Name = "Internal Energy";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = 1;
		break;
	case IO_RHO:
		Block.Label = "RHO ";
		Block.Name = "Density";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving =
		    (Comv2phys.Mass / pow(Comv2phys.Length, 3));
		break;
	case IO_HSML:
		Block.Label = "HSML";
		Block.Name = "Smoothing Length";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Length;
		break;
	case IO_VOL:
		Block.Label = "VOL ";
		Block.Name = "AREPO Cell Volume";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = pow(Comv2phys.Length, 3);
		break;
#ifdef BFLD
	case IO_BFLD:
		Block.Label = "BFLD";
		Block.Name = "Magnetic Field";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 3;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = 1;
		break;
#endif
	case IO_MASS:
		Block.Label = "MASS";
		Block.Name = "Particle Mass";
		for (i = 0; i < N_part_types; i++)
			Block.Npart[i] = Header.Npart[i];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Mass;
		break;
#ifdef VTURB
	case IO_VELT:
		Block.Label = "VELT";
		Block.Name = "Local Turbulent Velocity";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Vel;
		break;
	case IO_TNGB:
		Block.Label = "TNGB";
		Block.Name = "True Number Of Neighbours";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = INT;
		Block.Rmv_comoving = 1;
		break;
#endif
		/*Add above, not below !! */
	case IO_LASTENTRY:
		Block.Label = "LAST";
		Block.Name = "";
		Block.Val_per_element = 0;
		Block.Data_type = 0;
		break;
	}

	Block.Bytes_per_element = Block.Data_type * Block.Val_per_element;

	for (i = Block.Ntot = 0; i < N_part_types; i++)
		Block.Ntot += Block.Npart[i];

	return;
}

/*Fill P and Gas with data buffer 'fp'. 
 * */
void empty_comm_buffer(enum iofields blocknr, void *fp, size_t Pindex,
        size_t fpIndex)
{
	int i;

	switch (blocknr) {
	case IO_POS:
		for (i = 0; i < 3; i++)
			P[Pindex].Pos[i] = ((float *) fp)[fpIndex + i] 
                * Block.Rmv_comoving;
		break;
	case IO_VEL:
		for (i = 0; i < 3; i++)
			P[Pindex].Vel[i] = ((float *) fp)[fpIndex + i] 
                * Block.Rmv_comoving;
		break;
	case IO_ID:
		P[Pindex].ID = ((unsigned int *) fp)[fpIndex];
		break;
	case IO_U:
		Gas[Pindex].U = ((float *) fp)[fpIndex];
		break;
	case IO_RHO:
		Gas[Pindex].Rho = ((float *) fp)[fpIndex] * Block.Rmv_comoving;
		break;
	case IO_VOL:		/* Store VOL in HSML for now */
		Snap.Have_Arepo = 1;
	case IO_HSML:
		Gas[Pindex].Hsml = ((float *) fp)[fpIndex] * Block.Rmv_comoving;
		break;
#ifdef BFLD
	case IO_BFLD:
		for (i = 0; i < 3; i++)
			Gas[Pindex].Bfld[i] = ((float *) fp)[fpIndex + i];
		break;
#endif
	case IO_MASS:
		P[Pindex].Mass =
		    ((float *) fp)[fpIndex] * Block.Rmv_comoving;
		break;
#ifdef VTURB
	case IO_VELT:
		Gas[Pindex].VTurb =
		    ((float *) fp)[fpIndex] * Block.Rmv_comoving;
		break;
	case IO_TNGB:
		Gas[Pindex].TNgb = ((int *) fp)[fpIndex];
		break;
#endif
	case IO_LASTENTRY:  /*Add above, not below !! */
		break;
	}

	return;
}

/* 
 * Common input functions
 * */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define REAL 1
#define STRING 2
#define PINT 3
#define MAXTAGS 300

/* 
 * Reads a number of tags from an ascii file
 * the comment sign is  %
 * */
void read_param_file(char *filename)
{

	FILE *fd, *fdout;
	char buf[MAXLINELENGTH], buf1[MAXLINELENGTH], buf2[MAXLINELENGTH],
	    buf3[2 * MAXLINELENGTH];
	int i, j, nt;
	int id[MAXTAGS];
	void *addr[MAXTAGS];
	char tag[MAXTAGS][50];

	my_assert(sizeof(long long) == 8,
		  "Type `long long' is not 64 bit on this platform.");
	my_assert(sizeof(int) == 4,
		  "Type `int' is not 32 bit on this platform");
	my_assert(sizeof(float) == 4,
		  "Type `float' is not 32 bit on this platform");
	my_assert(sizeof(double) == 8,
		  "Type `double' is not 64 bit on this platform");

/* read parameter file on process 0 */
	if (ThisTask.Rank == 0) {
		nt = 0;

		strcpy(tag[nt], "Center_X");
		addr[nt] = &Param.Center[0];
		id[nt++] = REAL;

		strcpy(tag[nt], "Center_Y");
		addr[nt] = &Param.Center[1];
		id[nt++] = REAL;

		strcpy(tag[nt], "Center_Z");
		addr[nt] = &Param.Center[2];
		id[nt++] = REAL;

		strcpy(tag[nt], "Use_Barycenter");
		addr[nt] = &Param.Flag_Barycenter;
		id[nt++] = PINT;

		strcpy(tag[nt], "GridSize");
		addr[nt] = &Param.GridSize;
		id[nt++] = REAL;

		strcpy(tag[nt], "GridPoints");
		addr[nt] = &Param.NGrid;
		id[nt++] = PINT;

		strcpy(tag[nt], "Output_File");
		addr[nt] = &Param.Output_File;
		id[nt++] = STRING;

		strcpy(tag[nt], "Input_File");
		addr[nt] = &Param.Input_File;
		id[nt++] = STRING;

		strcpy(tag[nt], "N_IOTasks");
		addr[nt] = &Param.N_IOTasks;
		id[nt++] = PINT;

		strcpy(tag[nt], "Cosmology");
		addr[nt] = &Param.Cosmology;
		id[nt++] = PINT;

		strcpy(tag[nt], "NoClobber");
		addr[nt] = &Param.NoClobber;
		id[nt++] = PINT;

		strcpy(tag[nt], "UnitLength_in_cm");
		addr[nt] = &Unit.Length;
		id[nt++] = REAL;

		strcpy(tag[nt], "UnitMass_in_g");
		addr[nt] = &Unit.Mass;
		id[nt++] = REAL;

		strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
		addr[nt] = &Unit.Vel;
		id[nt++] = REAL;

#ifdef POWERSPECTRUM
		strcpy(tag[nt], "Nbins");
		addr[nt] = &Param.Pk_nbins;
		id[nt++] = PINT;

#endif
		if ((fd = fopen(filename, "r"))) {
			sprintf(buf, "%s%s", filename, "-usedvalues");
			if (!(fdout = fopen(buf, "w"))) {
				fprintf(stderr,
					"error opening file '%s' \n", buf);
				exit(1);
			} else {
				printf("\nReading Parameter file <%s>\n",
					filename);
				while (fgets(buf, MAXLINELENGTH, fd)) {

					if (sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
						continue;

					if (buf1[0] == '%')
						continue;

					for (i = 0, j = -1; i < nt; i++)
						if (strcmp(buf1, tag[i]) == 0) {
							j = i;
							tag[i][0] = 0;
							break;
						}

					if (j >= 0) {
						switch (id[j]) {
						case REAL:
							*((double *) addr[j]) = atof(buf2);
							fprintf(fdout, "%-35s%g\n"
                                    , buf1,*((double*)addr[j]));
							break;
						case STRING:
							strcpy((char *) addr[j], buf2);
							fprintf(fdout, "%-35s%s\n", buf1, buf2);
							break;
						case PINT:
							*((int *) addr[j]) = atoi(buf2);
                            fprintf(fdout, "%-35s%d\n", buf1, 
                                    *((int *) addr[j]));
							break;
						}
					}
				}
				fclose(fd);
				fclose(fdout);
				printf("\n");
			}
		} else 
			my_assert(0, "Parameter file not found");

		for (i = 0; i < nt; i++) {
			if (*tag[i]) {
				fprintf(stderr,
					"Value for tag '%s' missing in parameter file '%s'.\n",
					tag[i], filename);
				my_assert(0, "Parameter file incomplete");
			}
		}
	}

	MPI_Bcast(&Param, sizeof(Param), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Unit, sizeof(Unit), MPI_BYTE, 0, MPI_COMM_WORLD);

	return;
}

#undef REAL
#undef STRING
#undef INT
#undef MAXTAGS
