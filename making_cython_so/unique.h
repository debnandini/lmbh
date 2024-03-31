int main();
int *int_vector(int N);
void free_int_vector(int *v);
double **double_matrix(int N, int M);
void free_double_matrix(double **m, int N);
int **int_matrix(int N, int M);
void free_int_matrix(int **m, int N);
double *double_vector(int N);
void free_double_vector(double *v);
double ***double_tensor(int N, int M, int L);
void free_double_tensor(double ***t, int N, int M);

double logL(int M, int Nseg, double *cc, double *ss, double **CT, double **ST, double *SMT);

#define TPI 6.2831853071795865
#define year 3.15581498e7
