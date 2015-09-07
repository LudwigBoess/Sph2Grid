/* Common Function prototypes */
void *my_malloc(size_t size);
void my_free(void *ptr);
void *my_realloc(void *ptr, size_t size);
void Reallocate_P(long long, long long *, int);
void my_assert_full(int, const char *, const char *, const char *, int);

void print_compile_time_settings();

void null_final_operations();
void divide_npart_final_operations();

void center();
void domain_decomposition();
void set_workmode();

void init_grid();
void fill_grid();
void distribute_grid();
void test_grid();

void fft_grid();
void make_kGrid();
void powerspectrum();

void write_output();

/* Convert grid indices 3D to 1D */
size_t Idx(ptrdiff_t, ptrdiff_t, ptrdiff_t);
