
/* We set the workmode at compile time */
#ifdef SPH			/* e.g. Springel 2004 */
#define WORKMODE "SPH"
#define WMODE_INIT set_sph_factors()
#define WMODE_PART_SIZE sph_part_size
#define WMODE_EDGE_IDX sph_edge_idx
#define WMODE_WEIGHTS sph_weights
#define WMODE_FINAL_OPS sph_final_operations
#ifdef OVERSAMPLING
#define WMODE_MSG "Using 27 x oversampling of kernel"
#else
#define WMODE_MSG " "
#endif
#endif				/* SPH */

#ifdef NGP			/* Hockney & Eastwood */
#define WORKMODE "NGP"
#define WMODE_INIT
#define WMODE_PART_SIZE ngp_part_size
#define WMODE_EDGE_IDX ngp_edge_idx
#define WMODE_WEIGHTS ngp_weights
#define WMODE_FINAL_OPS divide_npart_final_operations
#define WMODE_MSG " "
#endif				/* NGP */

#ifdef CIC			/* Hockney & Eastwood */
#define WORKMODE "CIC"
#define WMODE_INIT
#define WMODE_PART_SIZE cic_part_size
#define WMODE_EDGE_IDX cic_edge_idx
#define WMODE_WEIGHTS cic_weights
#define WMODE_FINAL_OPS null_final_operations
#define WMODE_MSG " "
#endif				/* CIC */

#ifdef TSC			/* Hockney & Eastwood */
#define WORKMODE "TSC"
#define WMODE_INIT
#define WMODE_PART_SIZE tsc_part_size
#define WMODE_EDGE_IDX tsc_edge_idx
#define WMODE_WEIGHTS tsc_weights
#define WMODE_FINAL_OPS divide_npart_final_operations
#define WMODE_MSG " "
#endif				/* TSC */

#ifdef D20			/* Cui 2009 */
#define WORKMODE "D20 Wavelet"
#define WMODE_INIT set_dwave_vars()
#define WMODE_PART_SIZE dwave_part_size
#define WMODE_EDGE_IDX dwave_edge_idx
#define WMODE_WEIGHTS dwave_weights
#define WMODE_FINAL_OPS null_final_operations
#define WMODE_MSG  " "
#endif				/* D20 */

void divide_npart_final_operations();

/* SPH scheme */
double sph_part_size(long long);
void sph_edge_idx(double, double, double, double,
		  struct Cell_boundary_indices *);
double sph_weights(long long, double, double, double,
		   double, double, double);
void sph_final_operations();
void set_sph_factors();

/* Nearest Grid Point scheme */
double ngp_part_size(long long);
void ngp_edge_idx(double, double, double, double,
		  struct Cell_boundary_indices *);
double ngp_weights(long long, double, double, double,
		   double, double, double);

/* Cloud in Cell scheme */
double cic_part_size(long long);
void cic_edge_idx(double, double, double, double,
		  struct Cell_boundary_indices *);
double cic_weights(long long, double, double, double,
		   double, double, double);

/* Triangular Shaped Cloud scheme */
double tsc_part_size(long long);
void tsc_edge_idx(double, double, double, double,
		  struct Cell_boundary_indices *);
double tsc_weights(long long, double, double, double,
		   double, double, double);

/* Deaubechies D20 scheme */
void set_dwave_vars();
double dwave_part_size(long long);
void dwave_edge_idx(double, double, double, double
		    , struct Cell_boundary_indices *);
double dwave_weights(long long, double, double, double,
		     double, double, double);
