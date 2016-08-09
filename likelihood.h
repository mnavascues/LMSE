double p_delta_k(int delta, int k);
double p_DELTA_K(int DELTA[], int K[]);
double p_DELTA_j(int DELTA[], int j);
double p_j_theta_theta_tau(int j, double theta0, double theta1, double tau);
double p_j_theta(int j, double theta);
double p_DELTA_theta_theta_tau(int DELTA[],double theta0,double theta1,double tau);

double p_sample(double *punto);
double p_sample_no_H(double *punto);

void combination(int DELTA[],int maxK[],int K[],int counter,int j,double *pointer_sum);

double factorial(int M);

int max_mut_SMM(int pair);



extern int number_of_pop;
extern int total_sample_size;
extern int number_of_loci;
extern int **data;
extern double **delta_vector;
extern int *pop_size;
extern char **pop_name;
extern double *p;
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
extern int problem;
