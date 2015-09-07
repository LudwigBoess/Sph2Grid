/*Code Units*/
struct units {			/* For unit Conversion. */
	double Length;
	double Mass;
	double Vel;
	double Time;
	double Energy;
} Unit, Comv2phys;

double T(long long);
double Rho(long long);
double NumRho(long long);
double NumRhoProton(long long);	/* = n_th */

void set_units();
