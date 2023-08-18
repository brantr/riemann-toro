#include<stdio.h>
#include"routines.h"
#include"riemann_toro.h"
int main(int argc, char **argv)
{
	RiemannProblem RP;

	double stl[5];
	double str[5];
	double *sta;


    //toro p 129 test 1
 
    double rho_l = 1.0;
    double u_l   = 0.0;
    double p_l   = 1.0;
        //
    double rho_r  = 0.125;
    double u_r    = 0.0;
    double p_r    = 0.1;
    double gamma = 1.4;


	//density, velocity, pressure, gamma
	stl[0] = rho_l;
	stl[1] = u_l;
	stl[2] = p_l;
	stl[3] = gamma;
	stl[4] = gamma;

	str[0] = rho_r;
	str[1] = u_r;
	str[2] = p_r;
	str[3] = gamma;
	str[4] = gamma;

	sta = calloc_double_array(5);

	//RP.Test_Problem_Toro_1();
	//RP.Test_Problem_Toro_5();


	//RP.riemann(0,stl,str,sta);

	//printf("sta %f %f %f %f %f\n",sta[0],sta[1],sta[2],sta[3],sta[4]);



	FILE *fp;
	fp = fopen("toro_problem_1.txt","w");
	double t_final = 0.25;
	double s_min =  -2.0;
	double s_max =   2.0;
	double s;

	int i;
	int n_sample = 1000;

	for(i=0;i<n_sample;i++)
	{
		s = (s_max-s_min)*((double) i)/((double) (n_sample-1)) + s_min;

		RP.riemann(s,stl,str,sta);

		fprintf(fp,"%e\t%e\t%e\t%e\n",s*t_final+0.5,sta[0],sta[1],sta[2]);
	}
	fclose(fp);


	free(sta);

	return 0;
}
