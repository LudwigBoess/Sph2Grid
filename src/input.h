/*	Auxilliary Simulation properties  */

extern void read_param_file();

/* This Function fills the Particle structure
 * P and G, and sets in the Comv structure.
 * It should do so in parallel
 * */
extern void read_snapshot(char *);

/*Additional Stuff for GADGET-2 Input */

int swap;

enum iofields {
	IO_POS,
	IO_VEL,
	IO_ID,
	IO_U,
	IO_RHO,
	IO_HSML,
	IO_VOL,
	IO_MASS,
#ifdef BFLD
	IO_BFLD,
#endif
#ifdef VTURB
	IO_VELT,
	IO_TNGB,
#endif
#ifdef MACH
	IO_MACH,
#endif
	IO_LASTENTRY		/* Keep this entry at the end for termination */
};

void swap_Nbyte(char *, int, int);
size_t my_fread(void *, size_t, size_t, FILE *);
int find_block(FILE *, char *);
int find_files(char *);
void read_file(char *, int, int);
int read_gadget_block(void *, char *, FILE *, size_t);
int read_gadget_head(FILE *);
void set_block_prop(enum iofields);
void empty_comm_buffer(enum iofields, void *, size_t, size_t);

/*Input Header*/
struct Gadget_head {
	long long Npart[N_part_types];
	double Mass[N_part_types];
	double Time;
	double Redshift;
	long long Nall[N_part_types];
	long NumFiles;
	double BoxSize;
	double Omega0;
	double OmegaLambda;
	double HubbleParam;
	int FlagSfr;
	int FlagFeedback;
	int FlagCooling;
	int FlagAge;
	int FlagMetals;
} Header;

/* Block specifications */
struct blockdef {
	char *Label;
	char *Name;
	void *DataPtr;
	long long Npart[N_part_types];
	long long Ntot;
	int Val_per_element;
	size_t Data_type;
	size_t Bytes_per_element;
	double Rmv_comoving;
} Block;
