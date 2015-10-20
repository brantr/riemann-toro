/*! \file routines.c
 *  \brief Function definitions for utility routines.
 */
#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_spline.h>
#include<gsl/gsl_statistics.h>
#include<time.h>
#include"routines.h"


/*! \fn double double_log10_index(int i, int n, double xmin, double xmax)
 *  \brief Provides the i^th out of n log10 incremented value between xmin and xmax.
 *	   Useful for creating a ordinate array for an interpolation.
 */
double double_log10_index(int i, int n, double xmin, double xmax)
{
	return pow(10.0, (log10(xmax) - log10(xmin))*((double) i)/((double) (n-1)) + log10(xmin));
}
/*! \fn double double double_linear_index(int i, int n, double xmin, double xmax)
 *  \brief Provides the i^th out of n linear incremented value between xmin and xmax.
 *	   Useful for creating a ordinate array for an interpolation.
 */
double double_linear_index(int i, int n, double xmin, double xmax)
{
	return (xmax - xmin)*((double) i)/((double) (n-1)) + xmin;
}
/*! \fn void check_args(int argc, char **argv, int num_args)
 *  \brief Ensures that the number of command line arguements equals num_args.
 */
void check_args(int argc, char **argv, int num_args)
{
	if(argc!=num_args)
	{
		printf("Error: argc == %d, expected num_args == %d.\n",argc,num_args);
		fflush(stdout);
		exit(-1);
	}
}

/*! \fn int check_file(char fname[])
 *  \brief Check if a file exists.
 */
int check_file(char fname[])
{
        FILE *fp;
        int flag = 0;

        if(fp = fopen(fname,"r"))
        {
                flag = 1;
                fclose(fp);
        }
        return flag;
}


/*! \fn FILE *fopen_brant(char fname[], const char *mode)
 *  \brief Safe method for opening a FILE pointer.
 */
FILE *fopen_brant(char fname[], const char *mode)
{
	FILE *fp;

	if(!(fp = fopen(fname,mode)))
	{
		printf("Error opening %s.\n",fname);
		fflush(stdout);
		exit(-1);
	}

	return fp;
}
/*! \fn double *calloc_double_array(int n)
 *  \brief Safe method for callocing a double array
 */
double *calloc_double_array(int n)
{
	double *f;

	if(!(f = (double *) calloc(n,sizeof(double))))
	{
		printf("Error allocating array of size %d (%d).\n",n,(int) (n*sizeof(double)/sizeof(int)));
		fflush(stdout);
		exit(-1);
	}

	return f;
}


/*! \fn float *calloc_float_array(int n)
 *  \brief Safe method for callocing a float array
 */
float *calloc_float_array(int n)
{
	float *f;

	if(!(f = (float *) calloc(n,sizeof(float))))
	{
		printf("Error allocating array of size %d (%d).\n",n,(int) (n*sizeof(float)/sizeof(int)));
		fflush(stdout);
		exit(-1);
	}

	return f;
}


/*! \fn int *calloc_int_array(int n)
 *  \brief Safe method for callocing an int array
 */
int *calloc_int_array(int n)
{
	int *f;

	if(!(f = (int *) calloc(n,sizeof(int))))
	{
		printf("Error allocating array of size %d (%d).\n",n,(int) (n*sizeof(int)/sizeof(int)));
		fflush(stdout);
		exit(-1);
	}

	return f;
}

/*! \fn size_t *calloc_size_t_array(int n)
 *  \brief Safe method for callocing a size_t array
 */
size_t *calloc_size_t_array(int n)
{
	size_t *f;

	if(!(f = (size_t *) calloc(n,sizeof(size_t))))
	{
		printf("Error allocating array of size %d (%d).\n",n,(int) (n*sizeof(size_t)/sizeof(int)));
		fflush(stdout);
		exit(-1);
	}

	return f;
}
/*! \fn double max_three(double a, double b, double c)
 *  \brief Returns the max of 3 numbers */
double max_three(double a, double b, double c)
{
	return GSL_MAX_DBL(GSL_MAX_DBL(a,b),c);
}
/*! \fn float **two_dimensional_float_array(int n, int l)
 *  \brief Allocate a two dimensional (n x l) float array
 */
float **two_dimensional_float_array(int n, int l)
{
	float **x;

	x = new float *[n];
	for(int i=0;i<n;i++)
		x[i] = new float[l];

	for(int i=0;i<n;i++)
		for(int j=0;j<l;j++)
			x[i][j] = 0.0;
	return x;
}

/*! \fn double **two_dimensional_array(int n, int l)
 *  \brief Allocate a two dimensional (n x l) array
 */
double **two_dimensional_array(int n, int l)
{
	double **x;

	x = new double *[n];
	for(int i=0;i<n;i++)
		x[i] = new double[l];

	for(int i=0;i<n;i++)
		for(int j=0;j<l;j++)
			x[i][j] = 0.0;
	return x;
}
/*! \fn void deallocate_two_dimensional_array(double **x, int n, int l)
 *  \brief De-allocate a two dimensional (n x l) array
 */
void deallocate_two_dimensional_array(double **x, int n, int l)
{
	for(int i=0;i<n;i++)
		delete[] x[i];
	delete x;
}
void deallocate_two_dimensional_array(float **x, int n, int l)
{
	for(int i=0;i<n;i++)
		delete[] x[i];
	delete x;
}
/*! \fn double ***three_dimensional_array(int n, int l, int m)
 *  \brief Allocate a three dimensional (n x l x m) array 
 */
double ***three_dimensional_array(int n, int l, int m)
{
	double ***x;

	x = new double **[n];
	for(int i=0;i<n;i++)
	{
		x[i] = new double *[l];
		for(int j=0;j<l;j++)
		{
			x[i][j] = new double [m];
		}
	}

	return x;
}
/*! \fn void deallocate_three_dimensional_array(double ***x, int n, int l, int m)
 *  \brief De-allocate a three dimensional (n x l x m) array
 */
void deallocate_three_dimensional_array(double ***x, int n, int l, int m)
{
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<l;j++)
			delete[] x[i][j];
		delete[] x[i];
	}
	delete x;
}
/*! \fn double ***four_dimensional_array(int n, int l, int m, int k)
 *  \brief Allocate a four dimensional (n x l x m x k) array 
 */
double ****four_dimensional_array(int n, int l, int m, int p)
{
	double ****x;

	x = new double ***[n];
	for(int i=0;i<n;i++)
	{
		x[i] = new double **[l];
		for(int j=0;j<l;j++)
		{
			x[i][j] = new double *[m];
			for(int k=0;k<p;k++)
			{
				x[i][j][k] = new double [p];
			}
		}
	}

	return x;
}
/*! \fn void deallocate_four_dimensional_array(double ****x, int n, int l, int m, int p)
 *  \brief De-allocate a four dimensional (n x l x m x p) array
 */
void deallocate_four_dimensional_array(double ****x, int n, int l, int m, int p)
{
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<l;j++)
		{
			for(int k=0;k<p;k++)
				delete[] x[i][j][k];
			delete[] x[i][j];
		}
		delete[] x[i];
	}
	delete x;
}


/*! \fn int ***three_dimensional_int_array(int n, int l, int m)
 *  \brief Allocate a three dimensional (n x l x m) int array
 */
int ***three_dimensional_int_array(int n, int l, int m)
{
	int ***x;

	x = new int **[n];
	for(int i=0;i<n;i++)
	{
		x[i] = new int *[l];
		for(int j=0;j<l;j++)
		{
			x[i][j] = new int [m];
		}
	}

	return x;
}
/*! \fn void deallocate_three_int_dimensional_array(int ***x, int n, int l, int m)
 *  \brief De-allocate a three dimensional (n x l x m) int array.
 */
void deallocate_three_int_dimensional_array(int ***x, int n, int l, int m)
{
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<l;j++)
			delete[] x[i][j];
		delete[] x[i];
	}
	delete x;
}


/*! \fn int compare_doubles(const void *a, const void *b)
 *  \brief Function to compare doubles, for use with qsort.
 */
int compare_doubles(const void *a, const void *b)
{
	//compare doubles, from GNU C Library Ref. Manual
	const double *da = (const double *) a;
	const double *db = (const double *) b;
	return (*da > *db) - (*da < *db);
}

/*! \fn double time_in_seconds(clock_t A, clock_t B)
 *  \brief Returns time difference B-A in seconds for two clock_t
 */
double time_in_seconds(clock_t A, clock_t B)
{
	//stupid work around for icpc
	long bl = (long) B;
	long al = (long) A;
	double BA = (double) (bl-al);
	double time_check = BA;
	time_check = BA/CLOCKS_PER_SEC;
	return time_check;

	//return (double) (B-A)/CLOCKS_PER_SEC;
}





/*! \fn void create_log10_spline(double (*func)(double, void *), double *log10x, double *&log10y, int n, double *params, gsl_spline *&spline, gsl_interp_accel *&acc)
 *  \brief Routine to make a spline, interpolating in log10.
 */
void create_log10_spline(double (*func)(double, void *), double *log10x, double *&log10y, int n, double *params, gsl_spline *&spline, gsl_interp_accel *&acc)
{
	//log10x contains the locations at which
	//wish to evaluate the function func()

	//log10y is input unallocated
	//we will allocate it and it will
	//contain log10(func(x))

	//n is the size of the array of log10x and log10y
	
	//params is the array of parameters to pass to
	//func(x,params)

	//spline is the unallocated gsl_spline

	//acc is the unallocated gsl_interp_accel

	double xs;
	double ys;

	//allocate log10y
	log10y = calloc_double_array(n);

	//evaluate func at x locations
	for(int i=0;i<n;i++)
	{
		xs = pow(10.0,log10x[i]);

		ys = func(xs,params);

		if(ys<=0)
		{
			printf("func(%e) <= 0 (%e), cannot use log10 spline here.\n",xs,ys);
			fflush(stdout);
			exit(-1);
		}
		log10y[i] = log10(ys);
	}


	//allocate interpolants
	spline = gsl_spline_alloc(gsl_interp_cspline,n);

	acc    = gsl_interp_accel_alloc();

	//create interpoolation
	gsl_spline_init(spline, log10x, log10y, n);
}

/*! \fn void create_linear_spline(double (*func)(double, void *), double *x, double *&y, int n, double *params, gsl_spline *&spline, gsl_interp_accel *&acc)
 *  \brief Routine to create a spline, interpolated linear.
 */
void create_linear_spline(double (*func)(double, void *), double *x, double *&y, int n, double *params, gsl_spline *&spline, gsl_interp_accel *&acc)
{
	//x contains the locations at which
	//wish to evaluate the function func()

	//y is input unallocated
	//we will allocate it and it will
	//contain func(x)

	//n is the size of the array of x and y
	
	//params is the array of parameters to pass to
	//func(x,params)

	//spline is the unallocated gsl_spline

	//acc is the unallocated gsl_interp_accel

	//allocate y
	y = calloc_double_array(n);

	//evaluate func at x locations
	for(int i=0;i<n;i++)
		y[i] = func(x[i],params);


	//allocate interpolants
	spline = gsl_spline_alloc(gsl_interp_cspline,n);

	acc    = gsl_interp_accel_alloc();

	//create interpoolation
	gsl_spline_init(spline, x, y, n);
}

/*! \fn double array_max(double *x, int n)
 *  \brief Return the maximum of an array. */
double array_max(double *x, int n)
{
	/*double m = x[0];
	for(int i=1;i<n;i++)
	{
		m = GSL_MAX_DBL(m,x[i]);
	}*/
	double m = gsl_stats_max(x,1,n);
	return m;
}
/*! \fn double array_min(double *x, int n)
 *  \brief Return the minimum of an array. */
double array_min(double *x, int n)
{
	/*double m = x[0];
	for(int i=1;i<n;i++)
	{
		m = GSL_MIN_DBL(m,x[i]);
	}*/
	double m = gsl_stats_min(x,1,n);
	return m;
}

/*! \fn double vector_cross_product(double *x, double *y, int n); 
 *  \brief Find the cross product of x x y */
double *vector_cross_product(double *x, double *y, int ndim)
{
	double *cp;

	if(ndim==2)
	{
		cp = (double *) calloc(1,sizeof(double));
		cp[0] = x[0]*y[1] -  x[1]*y[0];
	}else{
		cp = (double *) calloc(3,sizeof(double));
		cp[0] = x[1]*y[2] -  x[2]*y[1];
		cp[1] = x[2]*y[0] -  x[0]*y[2];
		cp[2] = x[0]*y[1] -  x[1]*y[0];
	}


	return cp;
}

/*! \fn double vector_dot_product(double *x, double *y, int n); 
 *  \brief Find the dot product of x * y */
double vector_dot_product(double *x, double *y, int n)
{
	double dot = 0;

	for(int i=0;i<n;i++)
		dot+=x[i]*y[i];

	return dot;
}

/*! \fn double vector_magnitude(double *x, int n); 
 *  \brief Find the magnitude of x */
double vector_magnitude(double *x, int n)
{
	double dot = 0;

	for(int i=0;i<n;i++)
		dot+=x[i]*x[i];

	return sqrt(dot);
}

/*! \fn double **tensor_transformation(double **a, double **sigma, int ndim)
 *  \brief Apply transformation a to tensor sigma */
double **tensor_transformation(double **a, double **sigma, int ndim)
{
	

	double **result_A;
	double **result_B;
	double xt;

	result_A = two_dimensional_array(ndim, ndim);
	result_B = two_dimensional_array(ndim, ndim);

	//first do sigma a^T
	for(int i=0;i<ndim;i++)
	{
		for(int j=0;j<ndim;j++)
		{
			xt = 0.0;
			for(int l=0;l<ndim;l++)
				xt+=sigma[j][l]*a[l][j];	

			result_A[i][j] = xt;
		}
	}

	//now do a*(sigma a^T)
	for(int i=0;i<ndim;i++)
	{
		for(int j=0;j<ndim;j++)
		{
			xt = 0.0;
			for(int l=0;l<ndim;l++)
				xt+=a[j][l]*result_A[l][j];

			result_B[i][j] = xt;
		}
	}
	free(result_A);

	return result_B;
}
