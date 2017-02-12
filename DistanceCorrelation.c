#include<stdio.h>
#include<math.h>
#include<float.h>
#include<stdlib.h>
#include<time.h>

void dCOVi(double *x, double **df1,int irow, int jrow, int total,double ** dcormat1, int *byrow, int *dims,
              double *index, int *idx, double *DCOV);

double Akl(double **akl, double **A, int n);

/* functions in utilities.c */
 double **alloc_matrix(int r, int c);
 int    **alloc_int_matrix(int r, int c);
 void   free_matrix(double **matrix, int r, int c);
 void   free_int_matrix(int **matrix, int r, int c);
 void   permute(int *J, int n);
 void   roworder(double *x, int *byrow, int r, int c);
 void   Euclidean_distance(double *x, double **Dx, int n, int d);
 void   index_distance(double **Dx, int n, double index);
 void   vector2matrix(double *x, double **y, int N, int d, int isroworder);
 double ***alloc_matofmat(int s,int r, int c);
 void free_matofmat(double ***matrix, int s,int r, int c);


int main() {
int GENES = 12042;
int SAMPLES = 594;
double **df1;
double **dcormat1;
double x[2][5]= { {2,3,21,4,3 }, {12,13,121,14,13} };
double y[]= {12,23,1,45,32};
double dcov[] ={0,0,0,0};
int dims[]={SAMPLES,1,1,0,199};
int idx[]={7,5,1,2,6,4,3};
double pval = 1;
double index =1;
int byrow=1;
	
FILE *fp,*fp2,*fp3;
int i,j;
//allocating dataframe
df1 = (double**)malloc(GENES*sizeof(double*));
for( i =0;i<GENES;i++)
df1[i]=(double*)malloc(SAMPLES*sizeof(double));

//allocate dcormat1
dcormat1 = (double**)malloc(GENES*sizeof(double*));
for( i =0;i<GENES;i++)
dcormat1[i]=(double*)malloc(GENES*sizeof(double));

fp= fopen("C:/users/cumc/Desktop/tcgaov/fulldata/transposedfulldata.txt","r");
if (fp ==NULL) printf("open error");

for (i=0;i<GENES;i++)
 for (j=0;j<SAMPLES;j++)
  fscanf(fp, "%lf", &df1[i][j]);
i=0;
//start the process of computing distance matrices
dCOVi(df1[i],df1,i,i+1,GENES,dcormat1,&byrow, dims,&index, idx,dcov);
	
fp2=fopen("C:/users/cumc/Desktop/tcgaov/fulldata/dc.txt","w");
//printing the result file in a format understandable by R	
for( i=0;i<GENES;i++){
	if( i!=0)
	fprintf(fp2,"\n");
	for(j=0;j<GENES;j++){
fprintf(fp2,"%lf\t",dcormat1[i][j]); }}
fclose(fp2);
}

void dCOVi(double *x, double **df1,int irow, int jrow, int total,double ** dcormat1, int *byrow, int *dims,
              double *index, int *idx, double *DCOV) {
    /*  computes dCov(x,y), dCor(x,y), dVar(x), dVar(y)
        V-statistic is n*dCov^2 where n*dCov^2 --> Q
        dims[0] = n (sample size)
        dims[1] = p (dimension of X)
        dims[2] = q (dimension of Y)
        dims[3] = dst (logical, TRUE if x, y are distances)
        index : exponent for distance
        idx   : index vector, a permutation of sample indices
        DCOV  : vector [dCov, dCor, dVar(x), dVar(y)]
     */

    int    j, k, n, n2, p, q, dst;
    double ***Dx, ***A;
    double V;
    double *y;
    int ti,tj,i;

    n = dims[0];
    p = dims[1];
    q = dims[2];
    dst = dims[3];

    /* critical to pass correct flag dst from R */

    Dx = alloc_matofmat(total,n, n);
    printf("dx is succ allocated\n");
    A = alloc_matofmat(total,n,n);
    printf("A is succ allocated\n");
 
    for( i = 0;i<total;i++) {
       Euclidean_distance(df1[i], Dx[i], n, p);
       index_distance(Dx[i], n, *index);
       }
   printf("euc and index distances are calculated\n");

   for( i = 0;i<total;i++)
      Akl(Dx[i], A[i], n);
	
   printf("akl is succ calculated\n");
   printf("Time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);



  for(ti=0;ti<1;ti++)
  { printf("ti is %d\n",ti);
	  for(tj=ti+1;tj<total;tj++)
	{


    n2 = ((double) n) * n;

    /* compute dCov(x,y), dVar(x), dVar(y) */
    for (k=0; k<4; k++)
        DCOV[k] = 0.0;

    for (k=0; k<n; k++)
        for (j=0; j<n; j++) {
            DCOV[0] += A[ti][k][j]*A[tj][k][j];
            DCOV[2] += A[ti][k][j]*A[ti][k][j];
            DCOV[3] += A[tj][k][j]*A[tj][k][j];
        }

    for (k=0; k<4; k++) {
        DCOV[k] /= n2;
        if (DCOV[k] > 0)
            DCOV[k] = sqrt(DCOV[k]);
            else DCOV[k] = 0.0;
    }
    /* compute dCor(x, y) */
    V = DCOV[2]*DCOV[3];
    if (V > DBL_EPSILON)
        DCOV[1] = DCOV[0] / sqrt(V);
        else DCOV[1] = 0.0;

	}

}
printf("Time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
    free_matofmat(A,total, n, n);
    free_matofmat(Dx,total, n, n);

    return;
}

double Akl(double **akl, double **A, int n) {
    /* -computes the A_{kl} or B_{kl} distances from the
        distance matrix (a_{kl}) or (b_{kl}) for dCov, dCor, dVar
        dCov = mean(Akl*Bkl), dVar(X) = mean(Akl^2), etc.
    */
    int j, k;
    double *akbar;
    double abar;

    akbar = (double *)calloc(n, sizeof(double));
    abar = 0.0;
    for (k=0; k<n; k++) {
        akbar[k] = 0.0;
        for (j=0; j<n; j++) {
            akbar[k] += akl[k][j];
        }
        abar += akbar[k];
        akbar[k] /= (double) n;
    }
    abar /= (double) (n*n);

    for (k=0; k<n; k++)
        for (j=k; j<n; j++) {
            A[k][j] = akl[k][j] - akbar[k] - akbar[j] + abar;
            A[j][k] = A[k][j];
        }
    free(akbar);

    return(abar);
}


double **alloc_matrix(int r, int c)
{
    /* allocate a matrix with r rows and c columns */
    int i;
    double **matrix;
    matrix = (double **)calloc(r, sizeof(double *));
    for (i = 0; i < r; i++)
    matrix[i] =(double *) calloc(c, sizeof(double));
    return matrix;
}

double ***alloc_matofmat(int s,int r, int c)
{
    /* allocate a matrix with r rows and c columns */
    int i,j;
    double ***matrix;
    matrix = (double ***)malloc(s* sizeof(double **));

    for (i = 0; i < s; i++){

    matrix[i] =(double **) malloc(r* sizeof(double*));
	for(j=0;j<r;j++)
		{matrix[i][j] = (double *) malloc(c*sizeof(double));
		if(matrix[i][j]==NULL) printf("memory is NOT allocated\n");
		}
}


    return matrix;
}

int **alloc_int_matrix(int r, int c)
{
    /* allocate an integer matrix with r rows and c columns */
    int i;
    int **matrix;
    matrix = (int **)calloc(r, sizeof(int *));
    if(matrix == NULL)
    for (i = 0; i < r; i++)
    matrix[i] =(int *) calloc(c, sizeof(int));

    return matrix;
}

void free_matrix(double **matrix, int r, int c)
{
    /* free a matrix with r rows and c columns */
    int i;
    for (i = 0; i < r; i++) free(matrix[i]);
    free(matrix);
}

void free_matofmat(double ***matrix, int s,int r, int c)
{
    /* free a matrix with r rows and c columns */
    int i,j;
    for (i = 0; i < s; i++) {

		for( j =0;j<r;j++)
		free(matrix[i][j]);

		free(matrix[i]);
	}
    free(matrix);
}


void free_int_matrix(int **matrix, int r, int c)
{
    /* free an integer matrix with r rows and c columns */
    int i;
    for (i = 0; i < r; i++) free(matrix[i]);
    free(matrix);
}

void permute(int *J, int n)
{
    /*
       permute the first n integers of J
       if n is length(J), equivalent to R:
            J <- rev(sample(J, length(J), replace=FALSE))
    */
    int i, j, j0, m=n;
    for (i=0; i<n-1; i++) {
                  j = rand() % m;
		//j = floor(runif(0.0, (double) m));
        m--;
        j0 = J[j];
        J[j] = J[m];
        J[m] = j0;
    }
}

void vector2matrix(double *x, double **y, int N, int d, int isroworder) {
    /* copy a d-variate sample into a matrix, N samples in rows */
    int i, k;
    if (isroworder == 1) {
        for (k=0; k<d; k++)
            for (i=0; i<N; i++)
                y[i][k] = (*(x+i*d+k));
        }
    else {
        for (k=0; k<N; k++)
            for (i=0; i<d; i++)
                y[i][k] = (*(x+k*N+i));
        }
    return;
}


void roworder(double *x, int *byrow, int r, int c) {
    /*
      utility to convert a vector from column order to row order
      assume that x is r by c matrix as a vector in column order
    */
    int    i, j, k, n=r*c;
    double *y;
    if (*byrow == 1) return;
    y = (double *)calloc(n, sizeof(double));
    i = 0;
    for (j=0; j<r; j++) {
        for (k=0; k<n; k+=r) {
            y[i] = x[k+j];
            i++;
        }
    }
    for (i=0; i<n; i++)
        x[i] = y[i];
    free(y);
    *byrow = 1;
    return;
}


void distance(double **data, double **D, int N, int d) {
    /*
       compute the distance matrix of sample in N by d matrix data
       equivalent R code is:  D <- as.matrix(dist(data))
    */
    int    i, j, k;
    double dif;
    for (i=0; i<N; i++) {
        D[i][i] = 0.0;
        for (j=i+1; j<N; j++) {
            D[i][j] = 0.0;
            for (k=0; k<d; k++) {
                dif = data[i][k] - data[j][k];
                D[i][j] += dif*dif;
            }
            D[i][j] = sqrt(D[i][j]);
            D[j][i] = D[i][j];
        }
    }
    return;
}


void Euclidean_distance(double *x, double **Dx, int n, int d)
{
    /*
        interpret x as an n by d matrix, in row order (n vectors in R^d)
        compute the Euclidean distance matrix Dx
    */
    int i, j, k, p, q;
    double dsum, dif;
    for (i=1; i<n; i++) {
        Dx[i][i] = 0.0;
        p = i*d;
        for (j=0; j<i; j++) {
            dsum = 0.0;
            q = j*d;
            for (k=0; k<d; k++) {
                dif = *(x+p+k) - *(x+q+k);
                dsum += dif*dif;
            }
            Dx[i][j] = Dx[j][i] = sqrt(dsum);
        }
    }
}


void index_distance(double **Dx, int n, double index)
{
    /*
        Dx is an n by n Euclidean distance matrix
        if index NEQ 1, compute D^index
    */
    int i, j;

    if (fabs(index - 1) > DBL_EPSILON) {
        for (i=0; i<n; i++)
            for (j=i+1; j<n; j++) {
                Dx[i][j] = pow(Dx[i][j], index);
                Dx[j][i] = Dx[i][j];
            }
    }
}


void sumdist(double *x, int *byrow, int *nrow, int *ncol, double *lowersum)
{
    /*
       sum all pairwise distances between rows of x
       equivalent to this in R:  h <- sum(dist(x))
       x must be in row order: x=as.double(t(x))
    */

    int i, j, k, p, q, n=(*nrow), d=(*ncol);
    double sum, dsum, dif;
    if (*byrow == 0)
        roworder(x, byrow, n, d);
    sum = 0.0;
    for (i=1; i<n; i++) {
        p = i*d;
        for (j=0; j<i; j++) {
            dsum = 0.0;
            q = j*d;
            for (k=0; k<d; k++) {
                dif = *(x+p+k) - *(x+q+k);
                dsum += dif*dif;
            }
            sum += sqrt(dsum);
        }
    }
    (*lowersum) = sum;
}


