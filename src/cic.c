#include "allvars.h"

double cic_part_size(long long ipart)
{
	return (0.5 * Cell2kpc);
}

void
cic_edge_idx(double px, double py, double pz, double psize,
	     struct Cell_boundary_indices *idx)
{
	idx->Xmin = floor(px - psize);
	idx->Xmax = idx->Xmin + 2;

	idx->Ymin = floor(py - psize);
	idx->Ymax = idx->Ymin + 2;

	idx->Zmin = floor(pz - psize);
	idx->Zmax = idx->Zmin + 2;

	return;
}

double
cic_weights(long long ipart, double x, double y, double z,
	    double px, double py, double pz)
{
	return ((1 - fabs(x - px)) * (1 - fabs(y - py)) * 
			(1 - fabs(z - pz)));
}
