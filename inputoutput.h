int read_options(char *opt_file);
void read_input(char input_file[50]);
void read_head(FILE *fp);
void writes_head(FILE *foutput);
void read_1pop(FILE *fp);


extern int number_of_pop;
extern int total_sample_size;
extern int number_of_loci;
extern int **data;
extern double **delta_vector;
extern int *pop_size;
extern char **pop_name;
extern double *p;
extern double *mut;
extern double **mut95;
extern double q;
extern double mu;
extern int *num_of_pairs;
extern double **table_mutations;
extern int total_num_of_pairs;
extern int max_num_of_pairs;
extern int num_of_stored_pairs;
extern int pop;
extern int option;
extern FILE *flog;
//glabal variables to define options
extern int method_MCL;
extern int method_MM_c;
extern int method_MCL_c;
extern int method_MCLH;
extern int statistics;
extern double ref_for_stats[7];
extern char input_file[50];
extern char output_file[50];
extern int option_rates;
extern char input_rates[50];
extern int save_mismatch_option;
extern int bootstrap;

