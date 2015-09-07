/* Code Units */

#include "allvars.h"

void set_units();

void set_units()
{
	Unit.Time = Unit.Length / Unit.Vel;
	Unit.Energy = Unit.Mass * p2(Unit.Vel);

	if (!ThisTask.Rank)
		printf("Setting System of Units \n"
		       "   Unit Length = %g cm \n"
		       "   Unit Time   = %g sec\n"
		       "   Unit Mass   = %g g  \n\n", Unit.Length,
		       Unit.Time, Unit.Mass);

	return;
}

/*Defaults:  cm/kpc, g/(1e10 M_sol), cm/sec/(km/sec), */
struct units Comv2phys, Unit;

double Rho(long long ipart)
{
	return ((double) (Gas[ipart].Rho) * Unit.Mass /
		(Unit.Length * Unit.Length * Unit.Length));
}

double NumRho(long long ipart)
{
	return (Rho(ipart) * n2ne / (u_mol * m_p));
}

double NumRhoProton(long long ipart)
{
	return (Rho(ipart) * H_frac / m_p);
}

double T(long long ipart)
{				/*non-radiative for now */
	return (2.0 / 3.0 * Gas[ipart].U * Unit.Vel * Unit.Vel * m_p *
		mean_mol_weight / k_B);
}
