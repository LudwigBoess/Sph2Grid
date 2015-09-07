/* Debauchies Wavelet - (Cui et al 2008) */

#include "allvars.h"
#include "schemes.h"

#define KERNEL_TABLE_LENGTH 8192    // = 2^13 

#define ORDER  20			// filter 
#define FILTER_LENGTH (ORDER-1)
#define SPAN (ORDER*0.5)	// Daubechies Number = length of kernel
#define OFFSET (0)

void fill_kernel_table();
void calculate_LagDau_matrix();
double eval_daub_phi(long double);

static double Kernel_Table[KERNEL_TABLE_LENGTH];

/* Lagarias-Daubechies Matrices */
static long double  T0[FILTER_LENGTH][FILTER_LENGTH] = { { 0 } },
                    T1[FILTER_LENGTH][FILTER_LENGTH] = { { 0 } };


static const long double D20_filter[20] = { // Filter Coefficients 
	0.0377171575922413777705637759431856039930821429895071757833,
	0.2661221827938417898336528733330732036939096501263002545137,
	0.7455750714864667697776778760914390952756121571680221915727,
	0.9736281107336398913666196312898727233199645468633471108865,
	0.3976377417690173790221007363353688388876441712746318309064,
	-0.3533362017941126004987541055752144944722662932275311857642,
	-0.2771098787209663057125805876441331259190860176326387962663,
	0.1801274485333933326855370709079138216150003485139959905484,
	0.1316029871010700146473230122180511302559991975997229747905,
	-0.1009665711967794341394320358496627879130525910171935228534,
	-0.0416592480876016139620934376146842666922393506946278264366,
	0.0469698140973971216920437651915647222380073800011686176732,
	0.0051004369678144774653419022321615820492194836813768591199,
	-0.0151790023358564987876837378676726757263503242314241856677,
	0.001973325364963205188267088364746616645229218070440214677,
	0.0028176865901946762114567657265242814276812149045861789921,
	-0.0009699478398564109629691288141048182477091554832396375264,
	-0.0001647090060907779551031101592395127090433322701363872923,
	0.0001323543668511067663687628879693437522481530268050133757,
	-0.0000187584156275004083371169971592814738505996631128710397
};

static const long double D2_filter[4] = {
	0.6830127, 1.1830127,
	0.3169873, -0.1830127
};

double dwave_part_size(long long ipart)
{
	return (SPAN * Cell2kpc);
}

void dwave_edge_idx(double px, double py, double pz, double psize,
	struct Cell_boundary_indices *idx)
{

	idx->Xmin = floor(px) - ceil(OFFSET);
	idx->Xmax = floor(px) + SPAN - floor(OFFSET) + 1;

	idx->Ymin = floor(py) - ceil(OFFSET);
	idx->Ymax = floor(py) + SPAN - floor(OFFSET) + 1;

	idx->Zmin = floor(pz) - ceil(OFFSET);
	idx->Zmax = floor(pz) + SPAN - floor(OFFSET) + 1;

	return;
}

double dwave_weights(long long ipart, double x, double y, double z,
	double px, double py, double pz)
{
	size_t xidx, yidx, zidx;
	double dx, dy, dz;

	dx = (x - px + OFFSET);
	dy = (y - py + OFFSET);
	dz = (z - pz + OFFSET);

	if (dx >= SPAN || dx < 0 || dy >= SPAN || dy < 0 
	 || dz >= SPAN || dz < 0)
		return (0);

	xidx = floor(dx * KERNEL_TABLE_LENGTH / SPAN);
	yidx = floor(dy * KERNEL_TABLE_LENGTH / SPAN);
	zidx = floor(dz * KERNEL_TABLE_LENGTH / SPAN);

	return Kernel_Table[xidx] * Kernel_Table[yidx] * Kernel_Table[zidx];
}


void set_dwave_vars()
{
	calculate_LagDau_matrix();

	fill_kernel_table();

	return;
}

/* Lagarias Daubechies Cascading Algorithm 
 * Constructs wavelet to 64bit precision */
void fill_kernel_table()
{
	unsigned long long i;

	if (!ThisTask.Rank) 
        printf("Constructing Daubechies D20 Wavelet... "); fflush(stdout);

	long double dx = (long double) SPAN / KERNEL_TABLE_LENGTH;

#pragma omp parallel for private(i) shared(dx, Kernel_Table)
	for (i = 0; i < KERNEL_TABLE_LENGTH; i++)
		Kernel_Table[i] = eval_daub_phi(i * dx);

	if (!ThisTask.Rank)
		printf("done \n\n");

	return;
}

/* Daubechies Scaling function at x  
 * (Soman, Arathi 2009) */
double eval_daub_phi(long double x)
{
	const int nDigits = CHAR_BIT * sizeof(double) - 1;	// precision 

	int i, j, k;
	size_t xi;
	long double result[FILTER_LENGTH] = { 0 },
	    matrix[FILTER_LENGTH][FILTER_LENGTH] = { { 0 } }, 
	    lastMatrix[FILTER_LENGTH][FILTER_LENGTH] = { { 0 } };

	/* inverse (!) dyadic representation of decimal part of x */
	xi = floorl((x - floorl(x)) * pow(2L, nDigits + 1));

	for (i = 0; i < FILTER_LENGTH; i++) // start with a unity matrix 
		matrix[i][i] = 1L;

	for (int bit = 0; bit < nDigits; bit++) { // Cascade algorithm on bits
		for (i = 0; i < FILTER_LENGTH; i++)
			for (j = 0; j < FILTER_LENGTH; j++) {
				lastMatrix[i][j] = matrix[i][j];
				matrix[i][j] = 0L;
			}
		
        size_t bitmask = 1ULL << (nDigits - bit); // start at the end of xi 

		if (((xi & bitmask) >> (nDigits - bit)) & 0x01L) {

			for (i = 0; i < FILTER_LENGTH; i++) // xi[bit] = 1 
				for (j = 0; j < FILTER_LENGTH; j++)
					for (k = 0; k < FILTER_LENGTH; k++)
						matrix[i][j] += lastMatrix[i][k] * T1[k][j];
		} else {	

			for (i = 0; i < FILTER_LENGTH; i++) // xi[bit] = 0 
				for (j = 0; j < FILTER_LENGTH; j++)
					for (k = 0; k < FILTER_LENGTH; k++)
						matrix[i][j] += lastMatrix[i][k] * T0[k][j];
		}
	}
	
	for (i = 0; i < FILTER_LENGTH; i++) // row average 
		for (j = 0; j < FILTER_LENGTH; j++)
			result[i] += matrix[i][j];

	for (i = 0; i < FILTER_LENGTH; i++)
		result[i] /= (long double) FILTER_LENGTH;

	return result[(int) floorl(x)];
}

void calculate_LagDau_matrix()
{
	if (!ThisTask.Rank)
		printf("Constructing Lagarias Daubechies Matrices... ");

	for (int i = 0; i < FILTER_LENGTH; i++) {

		for (int j = 0; j < FILTER_LENGTH; j++) {

			int idx = 2 * i - j;
			
            if (idx <= FILTER_LENGTH && idx >= 0)
				T0[i][j] = D20_filter[idx];

			idx++;
			
            if (idx <= FILTER_LENGTH && idx >= 0)
				T1[i][j] = D20_filter[idx];
		}
	}

	if (!ThisTask.Rank)
		printf("done\n");

	return;
}
