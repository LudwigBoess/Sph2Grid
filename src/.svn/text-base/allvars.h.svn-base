/*Here we define common variables and functions, 
 * used in many routines PLEASE use sensible names. 
 * Example: m_p is better than prtn for proton mass */

#ifndef COMMON_H
#define COMMON_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <complex.h>		/* C99 complex types */

#include <mpi.h>

#include "config.h"

#ifdef OPENMP
#include <omp.h>
#endif

#ifdef FOURIERTRANSFORM
#include <fftw3.h>
#include <fftw3-mpi.h>
#endif

#include <gsl/gsl_math.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>

/* We need these virtually everywhere */
#include "cosmo.h"
#include "unit.h"
#include "proto.h"

/* Basic parameters */
#define N_part_types 6		/* Number of particle species */
#define MAXLINELENGTH 500	/* Max length of any I/O char buffer */

/*mathematical constants*/
#define pi 			M_PI
#define sqrt2		M_SQRT2
#define sqrt3       1.73205077648162841796875
#define fourpithirds 4.18879032135009765625

/*physical constants cgs*/
#define c			GSL_CONST_CGSM_SPEED_OF_LIGHT
#define e			GSL_CONST_CGSM_ELECTRON_CHARGE*GSL_CONST_CGSM_SPEED_OF_LIGHT
#define k_B 		GSL_CONST_CGSM_BOLTZMANN
#define m_p 		GSL_CONST_CGSM_MASS_PROTON
#define m_e			GSL_CONST_CGSM_MASS_ELECTRON

/* unit conversions */
#define barn2cgs	GSL_CONST_CGSM_BARN
#define eV2cgs		GSL_CONST_CGSM_ELECTRON_VOLT
#define GeV2cgs	    eV2cgs*1e-9
#define Msol2cgs    (double)(1.98892e33)
#define kpc2cgs 	(double)(3.08568025e21)
#define pc2cgs 	    (double)(3.08568025e18)
#define pix2cm		Param.XYSize/Param.XYPix*kpc2cgs
#define yr2sec      31556926

/* Chemistry */
#define H_frac		0.76	/* Hydrogen fraction */
#define He_frac	(1.0-H_frac)	/* Helium fraction */
#define u_mol		4.0/(5.0*H_frac+3.0)	/* Mean mol. weight in hydr. mass */
#define n2ne 		(H_frac+0.5*He_frac)/(2.0*H_frac+0.75*He_frac)
#define yHelium	He_frac / (4.0 *H_frac)
#define mean_mol_weight (1.0+4.0*yHelium)/(1.0+3.0*yHelium+1.0)

/* other */
#define Tcmb		(2.728*(1+Img.Redshift))	/* Temp of the CMB [K] */
#define Bcmb		(3.24516e-6*(1+Img.Redshift))	/* Magnetic field of the CMB [G] */

/* Macros */
#define my_assert(a,b) my_assert_full(a,b, __func__, __FILE__, __LINE__)

#define imin(a,b) ((a)>(b)?(b):(a))
#define imax(a,b) ((a)<(b)?(b):(a))

#define p2(a) ((a)*(a))
#define p3(a) ((a)*(a)*(a))

#define length3(a) sqrt(p2(a[0]) + p2(a[1]) + p2(a[2]))	/* this one is expensive ! */

#define Cell2kpc  (Param.GridSize / Param.NGrid)
#define Kpc2cell  (Param.NGrid / Param.GridSize)
#define Volume_cell (p3(Param.GridSize)/p3(Param.NGrid));
#define NGrid2 p2(Param.NGrid)
#define NGrid3 p3(Param.NGrid)

float fmaxf(float a, float b);
double fmax(double a, double b);
long double fmaxl(long double a, long double b);

float fminf(float a, float b);
double fmin(double a, double b);
long double fminl(long double a, long double b);

long labs(long i);
long long llabs(long long i);

double round(double x);

/* Profiling */
enum TimeMarks {
	CPU_READIN,
	CPU_DOMAINDECOMP,
	CPU_FILL,
    CPU_REDUCE,
	CPU_WORK,
	CPU_FFT,
	CPU_OUTPUT,
	CPU_LASTENTRY		/* Keep this entry at the end */
} mark;

struct Timings {
	double Zero;		/* Absolute Zeropoint */
	double Total[CPU_LASTENTRY];	/* Total Time for mark */
	double Start[CPU_LASTENTRY];	/* Start time for mark */
} Cpu;

void init_timing();
void finish_timing();
void start_timing(enum TimeMarks mark);
void stop_timing(enum TimeMarks mark);
double get_current_time();

/* Main Variables */
extern struct ParallelInfos {
	int Rank;		/* Rank of local Processor */
	int NTask;		/* Number of Processors */
	long long PartTotal;	/* Total Number of Particles on this Processor */
	long long Npart[N_part_types];	/* Npart stored locally */

	/* chunk of G struct stored on this CPU after grid reduction */
	size_t Nx;		/* stride of cut in k-space on this CPU */
	size_t X_start;		/* abs. start of cut in real-space on this CPU */
	size_t NDouble;		/* length of data in doubles */

#ifdef FOURIERTRANSFORM
	ptrdiff_t FFT_nx;
	ptrdiff_t FFT_x_start;
	ptrdiff_t FFT_nComplex;
#endif

} ThisTask;

#ifdef OPENMP
extern struct OpenMP_infos {
	int NThreads;		/* Number of openMP threads */
	int ThreadID;		/* Thread ID of this thread */
} Omp;
#pragma omp threadprivate(Omp)
#endif

extern struct InputParameters {	/*parameter from par file */
	double Center[3];
	double GridSize;	/* Extend of grid in 1D in [kpc] */
	size_t NGrid;		/* No of grid points along one dim. */
	int Flag_Barycenter;	/* Use Barycenter of Particles */
	int Cosmology;		/* Set in cosmo.c */
	int N_IOTasks;		/* Number of read tasks */
	int NoClobber;		/* Overwrite output file */
	char Output_File[MAXLINELENGTH];
	char Input_File[MAXLINELENGTH];
	int Workmode;		/* sets distribution scheme */
	size_t Pk_nbins;	/* for binned Powerspectrum */
} Param;

extern struct Snap_Properties {
	double Boxsize;
	double Redshift;
	double Time;
	double Masstab[N_part_types];
	long long PartTotal;
	long long Npart[N_part_types];
	int Have_Arepo;
} Snap;

/* Particle Data */
extern struct Particle_Data {
	float Pos[3];
	float Vel[3];
	unsigned long ID;
	int Type;
	float Mass;
} *P;

extern struct Gas_Data {
	float U;		/* Internal Energy */
	float Hsml;		/* Smoothing Length */
	float Rho;		/* Density */
#ifdef BFLD
	float Bfld[3];		/* Magnetic Field */
#endif
#ifdef VTURB
	float VTurb;		/* Turbulent Velocity */
	int TNgb;		/* True Number of Neighbours */
#endif
} *Gas;

/* Grid Variables & Handling */
extern struct Grid_Data {	/* Uses particle storage layout */
	double Mass;
	double Rho;
	double Npart;		/* Npart in cell */
#ifdef VEL
	double Vel[3];
#endif
#ifdef MOMENTUM
	double Mom[3];
#endif
#ifdef INTENERGY
	double U;
#endif
#ifdef BFLD
	double Bfld[3];		/* Magnetic Field */
#endif
#ifdef SCALAR_BFLD
	double B;		// Magnetic Field Strength
#endif
#ifdef VTURB
	double VTurb;
#endif
#ifdef DENSVEL
	double DensVel[3];	/* Density weighted velocity (Kritsuk+ 07) */
#endif
#ifdef VELDISPERSION
	double Sigma[3];
#endif
#ifdef SCALARVEL
	double VelScalar;
#endif
#ifdef VELFILTERED
	double VelFiltered[3];
#endif
} *Grid;

struct FT_Grid_Data {		/* Uses grid storage layout - no ifdefs */
	double complex *Mass;
	double complex *Rho;
	double complex *Vel[3];
	double complex *Mom[3];
	double complex *U;
	double complex *Bfld[3];
	double complex *B;
	double complex *VTurb;
	double complex *DensVel[3];
	double complex *VelScalar;
	double complex *VelFiltered[3];
	double *kVector[3];
} FTGrid;

struct powerspectra {
	double *Mass;
	double *Rho;
	double *Vel;
	double *Mom;
	double *U;
	double *Bfld;
	double *B;
	double *VTurb;
	double *DensVel;
	double *VelScalar;
	double *VelFiltered;
	double *k;
} Pk;

/* Needs to be global as dependent on WMODE */
struct Cell_boundary_indices {
	long Xmin, Xmax;
	long Ymin, Ymax;
	long Zmin, Zmax;
};

/* We use a global variable to hold the info
 * necessary to access the grid in an automated way.
 * It is set via the enum gridfields and 
 * set_gridfield_prop(). 
 * See init_grid()/setup.c for an example */
enum gridfields {		/* handle presence of grid fields */
	GRID_MASS,
	GRID_RHO,
	GRID_VEL,
	GRID_MOM,
	GRID_U,
	GRID_BFLD,
	GRID_B,
	GRID_VTURB,
	GRID_NPART,
	GRID_DENSVEL,
	GRID_SCALARVEL,
	GRID_VELDISPERSION,
	GRID_VELFILTERED,
	GRID_LASTENTRY		/* Add above */
};

/* we make pointer to pointers to be able to 
 * anonymously allocate arrays via something like
 * *GFld.Pkptr = my_malloc(nBytes) */
struct Gridfield_prop {		/* holds grid properties */
	char name[MAXLINELENGTH];
	double *ptr;		/* all grid values a dbl so far */
	double complex **FTptr;	/* points to FFT grids */
	double **Pkptr;		/* points to address holding Pk arr */
	int rank;		/* number of elements per field */
} GFld;

void set_gridfield_prop(enum gridfields);	/* sets GFld properties */

#endif
