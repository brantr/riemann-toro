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

	RP.Test_Problem_Toro_1();
	RP.Test_Problem_Toro_5();


	RP.riemann(0,stl,str,sta);

	free(sta);

	return 0;
}
