struct DATA{
  double
		*q,     // q[m] quality
    **s;    // s[m][t] price
};

struct PAR{
	// demand and price-equation parameters
  double
		beta,    // discount factor
		alpha,
    sig_alpha,
		*gamma,
    *alphai,
    phi,
		rho_p, // consumer expectation process
    sd_p,
		var_p;
};

struct RND{
	double
		**nu;			// nu[i][k]: unobserved heterogeneity for person i for char k
};

struct VAR{
	double
		**s,			// s[r][t]
    **c,      // c[m][t]
		***cp,			// cp[m][i][t]
		****vf;   // vf[m][i][t][r] value function for person i at time t & state r in market m
};

