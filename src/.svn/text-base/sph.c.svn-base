#include "allvars.h"
#include "schemes.h"

#define  NUMDIMS 3      /*!< For 3D-normalized kernel */
#define  KERNEL_COEFF_1  2.546479089470   
#define  KERNEL_COEFF_2  15.278874536822
#define  KERNEL_COEFF_3  45.836623610466
#define  KERNEL_COEFF_4  30.557749073644
#define  KERNEL_COEFF_5  5.092958178941
#define  KERNEL_COEFF_6  (-15.278874536822)
#define  NORM_COEFF      4.188790204786    

static double sph_kernel(double, double);
static double cell2kpc_sq;

double sph_part_size(long long ipart)
{
    return(Gas[ipart].Hsml);
}

void sph_edge_idx(double px, double py, double pz, double phsml
        , struct Cell_boundary_indices *idx)
{
    idx->Xmin = floor(px - phsml);
    idx->Xmax = ceil(px + phsml)+1;
    
    idx->Ymin = floor(py - phsml);
    idx->Ymax = ceil(py + phsml)+1;
    
    idx->Zmin = floor(pz - phsml) ;
    idx->Zmax = ceil(pz + phsml)+1;

    return;
}

double sph_weights(long long  ipart, double x, double y, double z, 
        double px, double py, double pz)
{
    	double weight = 0;

    	const double prho = Gas[ipart].Rho;
    	const double pmass = P[ipart].Mass;
    	const double phsml = sph_part_size(ipart);
    
#ifndef OVERSAMPLING 	// SPH definition 
    	double r2 = (x-px)*(x-px) + (y-py)*(y-py) + (z-pz)*(z-pz);

        weight = sph_kernel(phsml, r2); // r2 still in pixel units 

#else 	// 27 sampling points inside cell 
        const double x0 = x-px;
        
        const double y0 = y-py;

        const double z0 = z-pz;
    	
        for (int i=-1; i<2; i++)
    		for (int j=-1; j<2; j++)
    			for (int k=-1; k<2; k++){
				    double r2 = p2(x0+i*0.25) + p2(y0+j*0.25) + p2(z0+k*0.25);

                    weight += sph_kernel(phsml, r2);
			}

    	weight /= 27;
#endif	// OVERSAMPLING 

    	weight *= pmass / (prho);
    
    	return(weight);
}

#ifndef WC6
double sph_kernel(double hsml, double r2)
{

	double u;
	double hinv, hinv3;

    r2 *= cell2kpc_sq;

	if( r2 > hsml*hsml )
   	    return 0;

  	hinv = 1.0 / hsml;

  	hinv3 = hinv * hinv * hinv;
  
  	u = sqrt( r2 ) * hinv;

  	if(u < 0.5)
    	return(hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u));  
   	else 
    	return(hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u));
}
#else
double sph_kernel(double hsml, double r2)
{   
    r2 *= cell2kpc_sq;
	
	if( r2 > hsml*hsml )
   	    return 0;

	const double r = sqrt(r2);

	const double u= r/hsml;
    const double t = 1-u;

    return 1365.0/(64*pi)/p3(hsml) *t*t*t*t*t*t*t*t*(1+8*u + 25*u*u + 32*u*u*u);
}
#endif // WC6

void sph_final_operations()
{
    const ptrdiff_t nFloat_per_Cell = sizeof(*Grid)/sizeof(double);

    const double mean_npart_per_cell = (double)(Snap.PartTotal) / NGrid3;

    for (int block=0; block<GRID_LASTENTRY; block++) {

        set_gridfield_prop((enum gridfields)block);

        if (! GFld.rank)   // doesn't exist ? 
            continue;

        if (block == GRID_RHO) // doesn't need it
            continue;
        
        for(size_t i=0; i<p3(Param.NGrid); i++) // this is an approximation
            for(size_t j=0; j<GFld.rank; j++)
                GFld.ptr[i*nFloat_per_Cell+j] *= mean_npart_per_cell ; 
    }

    return;
}

void set_sph_factors()
{
    cell2kpc_sq = Cell2kpc*Cell2kpc;

    return;
}
