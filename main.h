
#include "inputoutput.h"
#include "Xatools.h"
#include "rates.h"
#include "rogers.h"
#include "likelihood.h"
#include "Hooke.h"

#define MAXSIZE 100
#define MAXLOCI 100

// GLOBAL VARIABLES
int number_of_pop;
int total_sample_size;
int number_of_loci;
int **data;
double **delta_vector;
int *pop_size;
char **pop_name;
double *p;
double q;
double mu;
int *num_of_pairs;
double **table_mutations;
int total_num_of_pairs;
int max_num_of_pairs;
int num_of_stored_pairs;
int pop;
int option=1;
FILE *flog;
int problem;
double **results;
double **table_stats;

//glabal variables to define options
int method_MCL;
int method_MM_c;
int method_MCL_c;
int method_MCLH;
int statistics;
double ref_for_stats[7];
char input_file[50]="data.txt";
char output_file[50]="results.csv";
char input_rates[50]="rates.txt";
int option_rates;
int save_mismatch_option;
int bootstrap;
double *mut;
double **mut95;

