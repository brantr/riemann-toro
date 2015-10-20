#ifndef   RIEMANN_TORO_BRANT
#define   RIEMANN_TORO_BRANT

//#define LOUD_RIEMANN

struct RiemannStructure
{
	double g1,g2,g3,g4,g5,g6,g7,g8;

	double rho_L;
	double u_L;
	double p_L;
	double c_sound_L;


	double rho_R;
	double u_R;
	double p_R;
	double c_sound_R;

	double gamma;

	double mpa;
};

class RiemannProblem
{
	public:
		RiemannStructure r;

		int riemann(double s, double stl[5], double str[5], double sta[5]);
		double guess_p_riemann(void);
		void   star_pu_riemann(double *pget, double *uget);
		void   prefun_riemann(double *f, double *fd, double p, double dk, double pk, double ck);
		void   sample_riemann(double s, double pm, double um, double *rho_s, double *u_s, double *p_s);

		void Test_Problem_Toro_1(void);
		void Test_Problem_Toro_2(void);
		void Test_Problem_Toro_3(void);
		void Test_Problem_Toro_4(void);
		void Test_Problem_Toro_5(void);
};





#endif  //RIEMANN_TORO_BRANT
