struct DATA {
	int
		***flag,    // indicate whether there's data
	    *i,   // individual info:  dummy org or ind (0 or 1)
	    *projwk;  // proj funding wk t end
	double
		***d,    // d[i][m][t] choice by person i at time t (continuous)
		***own_cumamt,     // own_cumamt[i][m][t] 
		***act_bound,     // own_cumamt[i][m][t] 
		***timeleft,     // timeleft[i][m][t] 
		****s,     // s[i][m][t][1...3] continuous state vector, 1: priorcumamt, 2 : priororgamt, 3: prior donor
		**q,      //proj information   1: goal for proj [m], 2: *logq,  3: *year, 4: *text , 5: * video, 6: * photo,7: * audio,
	//	****wkinfo,
		***vf;   // vf[m][i][t] given a specific observed state
};

struct PAR {
	// demand and price-equation parameters
	double
		//beta,     // discount factor
		*gamma,
		alpha,   // mean parameters
		sig_alpha,   // random-coefficient deviation
		*alphai, // individual-level parameter
		phi,      // decay parameter common across consumer & market
		rho_p,   // consumer expectation process
		sd_p, var_p;
};

struct RND{
	int
		i;
	double
		**nu,												// nu[i][k], unobserved heterogeneity for person i
		**s;
};

struct VAR {
	double
	//	**c,      // c[m][t] market specific value for person i at time t in market m
		****vf;   // vf[m][i][t][s] value function for person i at state s & time t in market m
};

struct LLH {
	double
		all,
		*i,
		**mi;  // ind[m][i] likelihood value for market m and person i
};

struct PSD {
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


