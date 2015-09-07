/* Test projected quantities */
#include "allvars.h"

void test_grid()
{
	size_t i;
    
    double sendRecvBuf[2] = { 0 };
	
    const size_t nGrid_local = ThisTask.NDouble;

    const size_t nPart_local = ThisTask.PartTotal;

	const double volume = p3(Param.GridSize / Param.NGrid);	// of one gridcell 

#ifdef MASS	 
	for (i = sendRecvBuf[0] = 0; i < nPart_local; i++)
		sendRecvBuf[0] += P[i].Mass;

	for (i = sendRecvBuf[1] = 0; i < nGrid_local; i++) // Grid already reduced
		sendRecvBuf[1] += Grid[i].Mass;

	MPI_Allreduce(MPI_IN_PLACE, sendRecvBuf, 2, MPI_DOUBLE, MPI_SUM,
		      MPI_COMM_WORLD);

	if (!ThisTask.Rank)
		printf("Total Mass in Particles <%g> \n"
		       "Total Mass on the Mass Grid  <%g> \n"
		       "Rel. Error <%g> \n\n", 
               sendRecvBuf[0] * Unit.Mass, sendRecvBuf[1] * Unit.Mass, 
               (sendRecvBuf[0] - sendRecvBuf[1]) / sendRecvBuf[0]);
#endif // MASS

#ifdef RHO	
   	for (i = sendRecvBuf[0] = 0; i < nPart_local; i++)
		sendRecvBuf[0] += P[i].Mass;

	for (i = sendRecvBuf[1] = 0; i < nGrid_local; i++) 
		sendRecvBuf[1] += Grid[i].Rho * volume;

	MPI_Allreduce(MPI_IN_PLACE, sendRecvBuf, 2, MPI_DOUBLE, MPI_SUM,
		      MPI_COMM_WORLD);

	if (!ThisTask.Rank)
		printf("Total Mass in Particles <%g> \n"
		       "Total Mass on the Density Grid  <%g> \n"
		       "Rel. Error <%g> \n\n", 
               sendRecvBuf[0] * Unit.Mass, sendRecvBuf[1] * Unit.Mass, 
               (sendRecvBuf[0] - sendRecvBuf[1]) / sendRecvBuf[0]);
#endif // RHO

#ifdef VEL
    for (i = sendRecvBuf[0] = 0; i < nPart_local; i++)
		sendRecvBuf[0] += p2(P[i].Vel[0]) + p2(P[i].Vel[1]) + p2(P[i].Vel[2]);

	for (i = sendRecvBuf[1] = 0; i < nGrid_local; i++) 
		sendRecvBuf[1] += p2(Grid[i].Vel[0]) + p2(Grid[i].Vel[1]) 
            + p2(Grid[i].Vel[2]);

	MPI_Allreduce(MPI_IN_PLACE, sendRecvBuf, 2, MPI_DOUBLE, MPI_SUM,
		      MPI_COMM_WORLD);

	if (!ThisTask.Rank)
		printf("Total v^2 energy in Particles <%g> \n"
		       "Total v^2 energy on the velocity Grid  <%g> \n"
		       "Rel. Error <%g> \n\n", 
               sendRecvBuf[0] * p2(Unit.Vel), sendRecvBuf[1] * p2(Unit.Vel), 
               (sendRecvBuf[0] - sendRecvBuf[1]) / sendRecvBuf[0]);
#endif // VEL

#ifdef MOMENTUM
    for (i = sendRecvBuf[0] = 0; i < nPart_local; i++)
		sendRecvBuf[0] += P[i].Mass * length3(P[i].Vel);

	for (i = sendRecvBuf[1] = 0; i < nGrid_local; i++) 
		sendRecvBuf[1] += length3(Grid[i].Mom);

	MPI_Allreduce(MPI_IN_PLACE, sendRecvBuf, 2, MPI_DOUBLE, MPI_SUM,
		      MPI_COMM_WORLD);

	if (!ThisTask.Rank)
		printf("Total momentum in Particles <%g> \n"
		       "Total momentum on the momentum Grid  <%g> \n"
		       "Rel. Error <%g> \n\n", 
               sendRecvBuf[0] * Unit.Mass * Unit.Vel, 
               sendRecvBuf[1] * Unit.Mass * Unit.Vel, 
               (sendRecvBuf[0] - sendRecvBuf[1]) / sendRecvBuf[0]);
#endif //MOMENTUM

#ifdef INTENERGY
    for (i = sendRecvBuf[0] = 0; i < nPart_local; i++)
		sendRecvBuf[0] += Gas[i].U;

	for (i = sendRecvBuf[1] = 0; i < nGrid_local; i++) 
		sendRecvBuf[1] += Grid[i].U;

	MPI_Allreduce(MPI_IN_PLACE, sendRecvBuf, 2, MPI_DOUBLE, MPI_SUM,
		      MPI_COMM_WORLD);

	if (!ThisTask.Rank)
		printf("Total internal Energy in Particles <%g> \n"
		       "Total internal Energy on the Grid  <%g> \n"
		       "Rel. Error <%g> \n\n", 
               sendRecvBuf[0], sendRecvBuf[1], 
               (sendRecvBuf[0] - sendRecvBuf[1]) / sendRecvBuf[0]);
#endif // INTENERGY

#ifdef BFLD
    for (i = sendRecvBuf[0] = 0; i < nPart_local; i++)
		sendRecvBuf[0] += p2(Gas[i].Bfld[0]) + p2(Gas[i].Bfld[1]) 
            + p2(Gas[i].Bfld[2]);
    sendRecvBuf[0] /= 8*pi;

	for (i = sendRecvBuf[1] = 0; i < nGrid_local; i++) 
		sendRecvBuf[1] += p2(Grid[i].Bfld[0]) + p2(Grid[i].Bfld[1]) 
            + p2(Grid[i].Bfld[2]);
    sendRecvBuf[1] /= 8*pi;

	MPI_Allreduce(MPI_IN_PLACE, sendRecvBuf, 2, MPI_DOUBLE, MPI_SUM,
		      MPI_COMM_WORLD);

	if (!ThisTask.Rank)
		printf("Total vector magnetic Energy in Particles <%g> \n"
		       "Total vector magnetic Energy on the Grid  <%g> \n"
		       "Rel. Error <%g> \n\n", 
               sendRecvBuf[0] * Unit.Energy, sendRecvBuf[1] * Unit.Energy, 
               (sendRecvBuf[0] - sendRecvBuf[1]) / sendRecvBuf[0]);
#endif // BFLD

#ifdef SCALAR_BFLD
	for (i = sendRecvBuf[0] = 0; i < nPart_local; i++)
		sendRecvBuf[0] += p2(Gas[i].Bfld[0]) + p2(Gas[i].Bfld[1]) 
            + p2(Gas[i].Bfld[2]);
    sendRecvBuf[0] /= 8*pi;

	for (i = sendRecvBuf[1] = 0; i < nGrid_local; i++) 
		sendRecvBuf[1] += p2(Grid[i].B);
    sendRecvBuf[1] /= 8*pi;

	MPI_Allreduce(MPI_IN_PLACE, sendRecvBuf, 2, MPI_DOUBLE, MPI_SUM,
		      MPI_COMM_WORLD);

	if (!ThisTask.Rank)
		printf("Total scalar magnetic Energy in Particles <%g> \n"
		       "Total scalar magnetic Energy on the Grid  <%g> \n"
		       "Rel. Error <%g> \n\n", 
               sendRecvBuf[0] * Unit.Energy, sendRecvBuf[1] * Unit.Energy, 
               (sendRecvBuf[0] - sendRecvBuf[1]) / sendRecvBuf[0]);
#endif // SCALAR_BFLD

#ifdef SCALARVEL
    for (i = sendRecvBuf[0] = 0; i < nPart_local; i++)
		sendRecvBuf[0] += length3(P[i].Vel);

	for (i = sendRecvBuf[1] = 0; i < nGrid_local; i++) 
		sendRecvBuf[1] += Grid[i].VelScalar;

	MPI_Allreduce(MPI_IN_PLACE, sendRecvBuf, 2, MPI_DOUBLE, MPI_SUM,
		      MPI_COMM_WORLD);

	if (!ThisTask.Rank)
		printf("Total vel power in Particles <%g> \n"
		       "Total scalar vel power on the Grid  <%g> \n"
		       "Rel. Error <%g> \n\n", 
               sendRecvBuf[0] * Unit.Vel, sendRecvBuf[1] * Unit.Vel, 
               (sendRecvBuf[0] - sendRecvBuf[1]) / sendRecvBuf[0]);
#endif //SCALARVEL

	MPI_Barrier(MPI_COMM_WORLD);

	return;
}
