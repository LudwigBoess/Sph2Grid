#include "allvars.h"

double ngp_part_size(size_t ipart)
{
	return Gas[ipart].Hsml;
}

void
ngp_edge_idx(double px, double py, double pz, double psize,
	     struct Cell_boundary_indices *idx)
{
	idx->Xmin = floor(px);
	idx->Xmax = idx->Xmin + 1;

	idx->Ymin = floor(py);
	idx->Ymax = idx->Ymin + 1;

	idx->Zmin = floor(pz);
	idx->Zmax = idx->Zmin + 1;

	return;
}

double ngp_weights(long long ipart, double x, double y, double z,
	    double px, double py, double pz)
{
	return 1;
}
