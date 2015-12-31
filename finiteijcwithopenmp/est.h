struct DATA{
	int
		***d;    // d[m][i][t] choice by person i at time t (0 or 1)
  double
		*q,      // q[m]
    **s,     // s[m][t] consumption normally distributed
		***vf;   // vf[m][i][t] given a specific observed state
};

struct PAR{
	// demand and price-equation parameters
  double
		beta,     // discount factor
		*gamma,
		alpha,   // mean parameters
		sig_alpha,   // random-coefficient deviation
		*alphai, // individual-level parameter
		phi,      // decay parameter common across consumer & market
		rho_p,   // consumer expectation process
		sd_p,var_p;
};

struct RND{
	int
		i;
	double
		**nu,												// nu[i][k], unobserved heterogeneity for person i
		**s;
};

struct VAR{
	double
    **c,      // c[m][t] market specific value for person i at time t in market m
		****vf;   // vf[m][i][t][s] value function for person i at state s & time t in market m
};

struct LLH{
  double
    all,
		*i,
    **mi;  // ind[m][i] likelihood value for market m and person i
};

struct PSD{
	int
		vmh,                       // current location in pvf array, cyclical from 1 to nvmh
		vmhmx,                     // if vmh has reached nvmh once, vmhmx=nvmh, otherwise vmhmx=vmh
		nvmh,                      // actual # pvf used to approximate pseudo-expected future values
    itmx,                      // max # of past pseudo-value functions to be used
    flgitmx,                   // flg for cycling included, = 0 no cycling
		ivmh;

	double
		**alphai, // thetai[i][1,...,1+n_mkt][vmh]
		**gamma,  // gamma[vmh]
		*phi,     // phi[vmh]
		**rs,     // rs[t][r]: history of state draw
		****pvf;	// pvf[vmh]: stored pseudo-value functions for buying decisions
		//*tkern,   // tkern[vmh]
		//*skern,   // skern[vmh]
		//mx_skern;
};


