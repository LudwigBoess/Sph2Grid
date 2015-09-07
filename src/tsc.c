#include "allvars.h"
#include "schemes.h"

double tsc_part_size(long long ipart)
{
	return (sqrt3 * Cell2kpc);
}

void
tsc_edge_idx(double px, double py, double pz, double psize,
	     struct Cell_boundary_indices *idx)
{
	idx->Xmin = floor(px) - 1;
	idx->Xmax = idx->Xmin + 3;

	idx->Ymin = floor(py) - 1;
	idx->Ymax = idx->Ymin + 3;

	idx->Zmin = floor(pz) - 1;
	idx->Zmax = idx->Zmin + 3;

	return;
}

double
tsc_weights(long long ipart, double x, double y, double z,
	    double px, double py, double pz)
{
	double weight = 1, dx, dy, dz;

	dx = fabs(px - x);
	dy = fabs(py - y);
	dz = fabs(pz - z);

	if (dx < 0.5) {
		weight *= 0.75 - dx * dx;
	} else {
		weight *= 0.5 * (1.5 - dx) * (1.5 - dx);
	}

	if (dy < 0.5) {
		weight *= 0.75 - dy * dy;
	} else {
		weight *= 0.5 * (1.5 - dy) * (1.5 - dy);
	}

	if (dz < 0.5) {
		weight *= 0.75 - dz * dz;
	} else {
		weight *= 0.5 * (1.5 - dz) * (1.5 - dz);
	}

	return (weight);
}
