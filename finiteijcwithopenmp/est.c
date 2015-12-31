
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <float.h>
#include <time.h>
//#include <sys/time.h>
#include "tools.h"
#include "est.h"
#include <omp.h>
#include "stdafx.h"

void pseudo_vf(const struct PAR *par, const struct DATA *data, struct PSD *psd, struct VAR *var, struct RND *rnd);
double kernel_app(const struct PAR *par, const struct DATA *data,const struct PSD *psd, const struct RND *rnd,
	const struct VAR *var, const int cm, const int ci, const int ct, const int flg);
void log_likelihood(struct LLH *llh, const struct PAR *par, const struct DATA *data, const struct PSD *psd,
	 struct VAR *var, const struct RND *rnd, const int flgi, const int flgm);
void true_vf(const struct PAR *par, const  struct DATA *data, const struct PSD *psd, struct VAR *var, const struct RND *rnd,
	const int cm, const int ci);
double exp_vf(const struct PAR *par, const struct DATA *data, const struct PSD *psd, const struct RND *rnd,
	const struct VAR *var, const int cm, const int ci, const int ct, const int flg);
void c_value(const struct PAR *par, const struct DATA *data, struct VAR *var,  int cm);

int
    n_param    =     8,  // # all parameters in param_input_sim.txt
	n_mkt      =    30,
	o_n_t      =    30,
	n_t        =   100,  // total # periods (terminal)
	n_person   =  100,
	n_theta    =     1,  // # hetero parameters
    n_gamma    =     2,  // # param for state transition density
	
	n_mh       =5000,
	n_pvmh     =  1000,
	n_rnd      =   100, // total # draws used for integration
	
	n_rj       =     0,

	bf_imh     =     0,

	flg_est;

double
	inv_2bwh;

const double
	euler_c    = 0.5772156649;


int main(void)
{
  //======================================================================
  // Declare variables
  //======================================================================
	int i,  m, t, d, k, r, imh, it, itemp1, itemp3; //itemp2, s, , flg, j,
	long seed=-123456;
	double dtemp1, dtemp2, dtemp3, dtemp4, dtemp5;
	time_t t1, t2, tst, tmt;
    inv_2bwh = 0.5 / 0.005;

	flg_est = 2; // =1 if full, =2 if pseudo

	printf("******************* OMP **********************\n");
	printf("%d\n", omp_get_num_threads());

	omp_set_num_threads(16);

	printf("******************* OMP **********************\n");
	printf("%d\n", omp_get_num_threads());

	double *param, maxs, mins, *jump, ln_prob_a, eps, acc, th_prob_a=1/3, *ap;

    struct DATA data;     // store data
    struct PAR par;       // store parameters
    struct PAR par_a;     // store accepted parameters
    struct PAR par_c;     // store candidate parameters
    struct LLH llh;       // store aggregate & individual llh
    struct LLH pllh;      // store aggregate & individual llh at previously accepted parameter
	struct VAR var;       // store variables
	struct RND rnd;       // store random numbers
	struct PSD psd;       // store pseudo-related variables (IJC, pseudo-FP)

    char param_name[8+1][30];

    FILE *inp, *outp;

  //======================================================================
  // Memory allocation
  //======================================================================
    param         = dvector(1, n_param);
	jump          = dvector(1, 4);
	ap            = dvector(1, 4);

	par.alphai    = dvector(1, n_person);
	par.gamma     = dvector(0, n_gamma);
	par_c.alphai  = dvector(1, n_person);
	par_c.gamma   = dvector(0, n_gamma);
	par_a.alphai  = dvector(1, n_person);
	par_a.gamma   = dvector(0, n_gamma);
	
	llh.mi        = dmatrix(1, n_mkt, 1, n_person);
    pllh.mi       = dmatrix(1, n_mkt, 1, n_person);
	llh.i         = dvector(1, n_person);
	pllh.i        = dvector(1, n_person);

	data.d        = iarray3d(1, n_mkt, 1, n_person, 1, n_t);
    data.q        = dvector(1, n_mkt);
    data.s        = dmatrix(1, n_mkt,  1, n_t);
	data.vf       = darray3d(1, n_mkt, 1, n_person, 1, n_t);

	rnd.nu        = dmatrix(1, n_person, 1, n_theta);
	rnd.s         = dmatrix(1, n_t, 1, n_rnd);

	psd.alphai    = dmatrix(1, n_person, 1, n_pvmh);
	psd.gamma     = dmatrix(0, n_gamma, 1, n_pvmh);
	psd.phi       = dvector(1, n_pvmh);
	psd.rs        = dmatrix(1, n_t, 1, n_pvmh);
	psd.pvf       = darray4d(1, n_mkt, 1, n_person, 1, n_t, 1, n_pvmh);
//	psd.tkern     = dvector(1, n_pvmh);
//	psd.skern     = dvector(1, n_pvmh);
  
    var.c         = dmatrix(1, n_mkt, 1, n_t);
	var.vf        = darray4d(1, n_mkt, 1, n_person, 2, n_t, 1, n_rnd);

  //============================================================================
  // Read in the starting parameter values
  //============================================================================
  printf("before reading starting parameter values.\n");

  inp = fopen("param_input.txt","r");
  i=1;

  while (i<=n_param){
    k = fscanf(inp, "%lf%s", &param[i], param_name[i]);
    printf("param[%d] = %f   %s.\n", i, param[i], param_name[i]);
    i++;
  }
	fclose(inp);

	par.beta      = par_c.beta      = par_a.beta      = param[1];
	par.gamma[0]  = par_c.gamma[0]  = par_a.gamma[0]  = param[2];
	par.gamma[1]  = par_c.gamma[1]  = par_a.gamma[1]  = param[3];
	par.alpha     = par_c.alpha     = par_a.alpha     = param[4];
	par.sig_alpha = par_c.sig_alpha = par_a.sig_alpha = param[5];
	par.phi       = par_c.phi       = par_a.phi       = param[6];
	par.rho_p     = par_c.rho_p     = par_a.rho_p     = param[7];
	par.sd_p      = par_c.sd_p      = par_a.sd_p      = param[8];

    printf("after reading starting parameter values.\n");

	par.var_p = par_c.var_p = par_a.var_p = par.sd_p * par.sd_p;

  //======================================================================
  // Read in data
  //====================================================================== 
  //--- Read data from data.txt ---------------------------
    printf("before reading data.\n");
  
	inp = fopen("data.txt", "r");  
	maxs=-DBL_MAX; mins= DBL_MAX;
	for (m=1;m<=n_mkt;m++){
		for (i=1;i<=n_person;i++){
			for (t=1;t<=o_n_t;t++){
        data.vf[m][i][t] = 0;
				data.s[m][t]     = 0;
        
				itemp1 = fscanf(inp, "%d%d%d%d%lf%lf%lf%lf%lf%lf%lf%lf",
					                    &m,&i,&t,&data.d[m][i][t],&data.s[m][t],&data.q[m],&data.vf[m][i][t],
															&dtemp1,&dtemp2,&dtemp3,&dtemp4,&dtemp5);

	      maxs = max(maxs, data.s[m][t]);
				mins = min(mins, data.s[m][t]);

				if (data.d[m][i][t]==1)
					break;

			} // t
		} // i
	} // m
	fclose(inp);

  // read nu.txt
  inp = fopen("nu.txt", "r");
	for (i=1;i<=n_person;i++){
		itemp3 = fscanf(inp, "%d%lf", &itemp1,&rnd.nu[i][1]);
	}
  fclose(inp);
  
  printf("after reading data.\n");

  //======================================================================
  // Preparation & Initialization
  //======================================================================

	// random variable generation for integration (for now use monte carlo)
//	dtemp1 = 0.2*abs(maxs-mins);
//	maxs += dtemp1;
//	mins -= dtemp1;
//	printf("maxs=%f, mins=%f\n",maxs,mins);
//	for (r=1;r<=n_rnd;r++)
//    for (t=1;t<=n_t;t++)
//			for (s=2;s<=n_state;s++){ // s=1 is state dependence
//	      rnd.s[r][t][s] = (maxs-mins)*ran1(&seed) + mins;
//				rnd.s[r][t][s] = 200*ran1(&seed)-100;
//			}

	inp = fopen("state.txt", "r");
	for (r=1;r<=n_rnd;r++){
		for (t=1;t<=n_t;t++){
			itemp3 = fscanf(inp, "%d%d%lf", &r, &t, &rnd.s[t][r]);
		}
	}
	fclose(inp);



  // initialization of pseudo-VF
	for (m=1;m<=n_mkt;m++)
		for (i=1;i<=n_person;i++)
			for (t=1;t<=n_t;t++)
				for (it=1;it<=n_pvmh;it++)
					psd.pvf[m][i][t][it] = 0;

	for (i=1;i<=n_person;i++){
		par.alphai[i] = par_c.alphai[i] = par_a.alphai[i] = par.alpha + par.sig_alpha * gasdev(&seed);
//		par.alphai[i] = par_c.alphai[i] = par_a.alphai[i] = par.alpha + par.sig_alpha * rnd.nu[i][1];
	}

	for (k=1;k<=4;k++){
		jump[k] = 0.01;
		ap[k]   = 0;
	}

  psd.vmh     = 1;       // location indicator of past pseudo-values array
  psd.vmhmx   = 1;       // max until mcmc iteration reaches n_pvmh
	psd.nvmh    = n_pvmh;  // set this to be < n_pvmh if using a subset of past pseudo-values to do approximation
	psd.itmx    = 0;       // max when involving cycling
	psd.flgitmx = 0;       // whether we need cycling
	psd.ivmh    = 1;       // pick one consumer each iteration to compute pseudo-value functions
	rnd.i       = 1;       // indicator for recycling random draws

	printf("after initialization\n");

  //======================================================================
  // MCMC Loop
  //======================================================================
  printf("starting MCMC loop\n\n");
  outp    = fopen("hist_param.txt", "w");
  t1 = gettimeofday_sec();
  tst = t1;

  for (imh=1; imh<=n_mh; imh++)
  {

		if (imh==1 || flg_est==2)
		{
		  // Log-likelihood functions at most recent H^r
			log_likelihood(&pllh, &par, &data, &psd, &var, &rnd, 0, 0);
		}
    
		// substituting the previously accepted values
		par.gamma[0]  = par_a.gamma[0];
		par.gamma[1]  = par_a.gamma[1];
		par.alpha     = par_a.alpha;
		par.sig_alpha = par_a.sig_alpha;
		for (i = 1; i <= n_person; i++)
		{
			par.alphai[i] = par_a.alphai[i];
		}
        par.phi = par_a.phi;

    
		if (flg_est == 2 && imh % 100 == 0)
		{
			dtemp1 = kernel_app(&par, &data, &psd, &rnd, &var, 1, 1, 2, 0);
			true_vf(&par, &data, &psd, &var, &rnd, 1, 1);
			dtemp2 = exp_vf(&par, &data, &psd, &rnd, &var, 1, 1, 2, 0);
			printf("t=1,vf=%f, tvf=%f, dvf=%f\n",dtemp1,dtemp2,data.vf[1][1][1]);
		}

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Update the pvf loop index given vmh,vmhmx,nvmh
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (flg_est == 2)
	{
      if (psd.vmhmx-1 < psd.nvmh && psd.vmhmx < n_pvmh)
	  { // case1 - vmhmx hasn't reached n_pvmh & vmh-1<nvmh
        psd.itmx    = psd.vmhmx - 1; // equivalent to vmh-1 (vmh has been incremented after storing most recent pseudo-value, so upto vmh-1)
        psd.flgitmx = 0;             // flg for cycling included, = 0 no cycling
      }else
	  { // case2&3, we can use nvmh pvfs, but may need cycling when indexing (case2 no cycling, case3 cycling)
        psd.itmx = psd.nvmh;
        if (psd.vmh-1 < psd.nvmh) psd.flgitmx = 1; // note vmh has been incremented after storing, so we don't use pvf at vmh
        else psd.flgitmx = 0;                      // if =1, we need to use vmh-1,...,1,n_pvmh,n_pvmh-1,....,n_pvmh-(nvmh-(vmh-1)+1)
      }
    }

if (imh>bf_imh){

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Draw population parameters (theta[1]-theta[1+n_mkt, sigma[1], sigma[2])
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	  // Draw mean parameter:
	  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		dtemp1 = 0;
		for (i=1;i<=n_person;i++){
			dtemp1 += par.alphai[i];
		}
		dtemp1 /= n_person;

		// draw from the posterior distribution (normal)
		par.alpha   = par_a.alpha = dtemp1 + par.sig_alpha/sqrt(n_person)*gasdev(&seed);
		par_c.alpha = par.alpha;

	  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	  // Draw sd parameter
	  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		dtemp1 = 0;
		for (i=1;i<=n_person;i++)
			dtemp1 += (par.alphai[i]-par.alpha)*(par.alphai[i]-par.alpha);
		dtemp1 /= n_person;

		// draw from the posterior distribution (Inverted gamma)
		dtemp2 = 0;
		for (d=1;d<=n_person+1;d++){
			dtemp3 = gasdev(&seed);
			dtemp2 += dtemp3*dtemp3;
		}
		dtemp2 /= (1+n_person*dtemp1);

		par.sig_alpha   = par_a.sig_alpha = sqrt(1/dtemp2);
		par_c.sig_alpha = par.sig_alpha;

		
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Draw individual parameters
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		pllh.all = 0;
  
		for (i=1;i<=n_person;i++){

  	  //------------------------------------------------------------------
  	  // Draw candidate parameter values, random-walk proposal distribution
  	  // & Compute log-likelihood at candidate parameter vector
  	  //------------------------------------------------------------------
      par.alphai[i] = par.alpha + par.sig_alpha * gasdev(&seed);
      par_c.alphai[i] = par.alphai[i];
      
  		log_likelihood(&llh, &par, &data, &psd, &var, &rnd, 0, i);

  	  //------------------------------------------------------------------
  	  // M-H acceptance or rejection
  	  //------------------------------------------------------------------
  	  ln_prob_a = min(llh.i[i]-pllh.i[i],0);
			//printf("llh.i[%d]=%f, pllh.i[%d]=%f\n",i,llh.i[i],i,pllh.i[i]);

  	  eps = ran1(&seed);
  	  if (log(eps)<=ln_prob_a){
  			par_a.alphai[i] = par.alphai[i];
				pllh.i[i] = llh.i[i];
        for (m=1;m<=n_mkt;m++){
          pllh.mi[m][i] = llh.mi[m][i];
				}
  	  }else{
  			par.alphai[i] = par_a.alphai[i];
  		}

			pllh.all += pllh.i[i]; // this is for the next hetero parameter
     
    } // i
	
	
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Draw homogeneous parameters (gamma)
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		for (k=0;k<=1;k++){

			//------------------------------------------------------------------
			// Draw candidate parameter values, random-walk proposal distribution
			// & Compute log-likelihood at candidate parameter vector
			//------------------------------------------------------------------
			par.gamma[k] = par.gamma[k] + jump[k+2] * gasdev(&seed);
			par_c.gamma[k] = par.gamma[k];

			// Compute log-likelihood at candidate parameter vector
			log_likelihood(&llh, &par, &data, &psd, &var, &rnd, 0, 0);
		
			//------------------------------------------------------------------
			// M-H acceptance or rejection
			//------------------------------------------------------------------
			ln_prob_a = min(llh.all-pllh.all,0);
		
			eps = ran1(&seed);
			if (log(eps)<=ln_prob_a){
				par_a.gamma[k] = par.gamma[k];
		    pllh.all = llh.all;
	      for (i=1;i<=n_person;i++)
	        pllh.i[i] = llh.i[i];
				ap[k+2]++;
			}else{
				par.gamma[k] = par_a.gamma[k];
			}

			// adjust the sd for random-walk proposal function based on acceptance probability
			acc = exp(ln_prob_a);
			if (acc>th_prob_a) jump[k+2] = 0.01;
			else jump[k+2] = 0.005;
		} // gamma k loop

		
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Draw homogeneous parameters (phi)
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		//------------------------------------------------------------------
		// Draw candidate parameter values, random-walk proposal distribution
		// & Compute log-likelihood at candidate parameter vector
		//------------------------------------------------------------------
		par.phi   = par.phi + jump[4] * gasdev(&seed);
    while (par.phi<0 || par.phi>1){
      par.phi   = par.phi + jump[4] * gasdev(&seed);
    }
		par_c.phi = par.phi;

		// Compute log-likelihood at candidate parameter vector
		log_likelihood(&llh, &par, &data, &psd, &var, &rnd, 0, 0);
		
		//------------------------------------------------------------------
		// M-H acceptance or rejection
		//------------------------------------------------------------------
		ln_prob_a = min(llh.all-pllh.all,0);
		//printf("llh=%f, pllh=%f\n",llh.all,pllh.all);
		
		eps = ran1(&seed);
		if (log(eps)<=ln_prob_a){
			par_a.phi = par.phi;
      pllh.all = llh.all;
      for (i=1;i<=n_person;i++)
        pllh.i[i] = llh.i[i];
			ap[4]++;
		}else{
			par.phi = par_a.phi;
		}

		// adjust the sd for random-walk proposal function based on acceptance probability
		acc = exp(ln_prob_a);
		if (acc>th_prob_a) jump[4] = 0.005;
		else jump[4] = 0.001;
	

}

	  //------------------------------------------------------------------
	  // Compute pseudo-value function at candidate parameter vector
	  //------------------------------------------------------------------
    if (flg_est==2)
      pseudo_vf(&par_c, &data, &psd, &var, &rnd);

    // write out the accepted parameter values in iteration [imh]
    if (imh % 1 == 0){
      tmt = gettimeofday_sec(); // take time here so it won't include true vf computation
      fprintf(outp, "%d\t%f",imh,pllh.all);
      fprintf(outp, "\t%f\t%f\t%f\t%f\t%f",par.gamma[0],par.gamma[1],par.alpha,par.sig_alpha,par.phi);
	  if (flg_est == 2) {
		  dtemp1 = kernel_app(&par, &data, &psd, &rnd, &var, 1, 1, 2, 0);
		  true_vf(&par, &data, &psd, &var, &rnd, 1, 1);
		  dtemp2 = exp_vf(&par, &data, &psd, &rnd, &var, 1, 1, 2, 0);
		  fprintf(outp, "\t%f\t%f", dtemp1, dtemp2);
	  }
      fprintf(outp, "\t%f\n", difftime(tmt,tst));
    }
    
    // print param
	if (imh % 1 == 0) {
		printf("\niteration   = %d / %d, flg_est = %d\n", imh, n_mh, flg_est);
		printf("llh	 	      = %.20f\n", pllh.all);

		printf("par_a.gamma[0]  =%f, par_c.gamma[0]  =%f,  ap = %f\n", par_a.gamma[0], par_c.gamma[0], ap[2] / imh);
		printf("par_a.gamma[1]  =%f, par_c.gamma[1]  =%f,  ap = %f\n", par_a.gamma[1], par_c.gamma[1], ap[3] / imh);
		printf("par_a.alpha     =%f, par_c.alpha    =%f, ap = %f\n", par_a.alpha, par_c.alpha, ap[1] / imh);
		printf("par_a.sig_alpha =%f, par_c.sig_alpha =%f\n", par_a.sig_alpha, par_c.sig_alpha);
		printf("par_a.phi       =%f, par_c.phi       =%f,  ap = %f\n", par_a.phi, par_c.phi, ap[4] / imh);
		printf("vmh = %d, vmhmx = %d, nvmh = %d\n", psd.vmh, psd.vmhmx, psd.nvmh);

		t2 = gettimeofday_sec();
		printf("sec (iter %d) = %f\n\n", imh, difftime(t2, t1));
		t1 = t2;

	}
    
  } // end of imh loop

	fclose(outp);

	printf("end!\n");

  //======================================================================
  // Free up Memory
  //======================================================================
    free_dvector(param,         1, n_param);
	free_dvector(jump,          1, 4);

	free_dvector(par.alphai,		1, n_person);
	free_dvector(par.gamma,     0, n_gamma);
	free_dvector(par_c.alphai,	1, n_person);
	free_dvector(par_c.gamma,   0, n_gamma);
	free_dvector(par_a.alphai,	1, n_person);
	free_dvector(par_a.gamma,   0, n_gamma);

	free_dmatrix(llh.mi,        1, n_mkt, 1, n_person);
    free_dmatrix(pllh.mi,       1, n_mkt, 1, n_person);
	free_dvector(llh.i,         1, n_person);
	free_dvector(pllh.i,        1, n_person);

	free_iarray3d(data.d,       1, n_mkt, 1, n_person, 1, n_t);
	free_dvector(data.q,        1, n_mkt);
    free_dmatrix(data.s,        1, n_mkt, 1, n_t);
    free_darray3d(data.vf,      1, n_mkt, 1, n_person, 1, n_t);

	free_dmatrix(rnd.nu,        1, n_person, 1, n_theta);
	free_dmatrix(rnd.s,         1, n_t, 1, n_rnd);

	free_dmatrix(psd.alphai,    1, n_person, 1, n_pvmh);
  	free_dmatrix(psd.gamma,     0, n_gamma, 1, n_pvmh);
    free_dvector(psd.phi,       1, n_pvmh);
	free_dmatrix(psd.rs,        1, n_t, 1, n_pvmh);
	free_darray4d(psd.pvf,      1, n_mkt, 1, n_person, 1, n_t, 1, n_pvmh);
	//free_dvector(psd.tkern,     1, n_pvmh);
	//free_dvector(psd.skern,     1, n_pvmh);

    free_dmatrix(var.c,         1, n_mkt, 1, n_t);
	free_darray4d(var.vf,       1, n_mkt, 1, n_person, 2, n_t, 1, n_rnd);
  
	return 0;
}


/*******************************************************************************
 Function: log_likelihood
 *******************************************************************************/
void log_likelihood(struct LLH *llh, const struct PAR *par, const struct DATA *data, const struct PSD *psd,
	 struct VAR *var, const struct RND *rnd, const int flgm, const int flgi)
{
  //======================================================================
  // Declare variables
  //======================================================================
	int m, i, t, ibegin, iend; // mbegin, mend;
  //int tid, nThreads;
  double av,u0,u1,pa1; //vf0,vf1,
  double **avf = dmatrix(1, n_mkt, 1, n_t);

//  printf("inside log_likelihood\n");
  
	//FILE *outp;

	if (flgi>0){
		ibegin = flgi;
		iend   = flgi;
	}else{
		ibegin = 1;
		iend = n_person;
	}

	for (int m=1;m<=n_mkt;m++){
		// c_value - compute a sequence of var->c[m][i][t=1,...,n_t]
		c_value(par, data, var, m);
	}


	//printf("Before setting the pragma, Number of threads: %d\n", omp_get_num_threads());
	#pragma omp parallel for schedule(dynamic) private(i, m, t,u0,u1,av,pa1)
	for ( i=ibegin;i<=iend;i++)
	{
		/*tid = omp_get_thread_num();
		printf("Running on thread %d\n", tid);

		if (tid == 0)
		{
			nThreads = omp_get_num_threads();
			printf("Total number of threads: %d\n", nThreads);
		}*/

		llh->i[i] = 0;

		for (m=1;m<=n_mkt;m++)
		{
		    llh->mi[m][i] = 0;

			// HETEROGENEOUS CASE
			// if using true value function
			if (flg_est == 1) 
			{
				true_vf(par, data, psd, var, rnd, m, i);
			}

			for (t = 1; t <= o_n_t; t++)
			{
				// HETEROGENEOUS CASE
				if (t == n_t)
				{
					av = 0;
				}
				else
				{
					if (flg_est == 1)
					{
						av = exp_vf(par, data, psd, rnd, var, m, i, t + 1, 0);
					}
					else if (flg_est == 2)
					{
						av = kernel_app(par, data, psd, rnd, var, m, i, t + 1, 0);
					}
				}
			//	
			//	if (isnan(av)) {
			//	//	int waitme;
			//		printf("\n out of kernal hit:  *av = %d,flgi=%d", av, flgi);
			////		scanf("input var ?%d", &waitme);
			//	}

				u0 = par->beta * av; // HETERO
//				u0 = par->beta * avf[m][t]; // HOMO

				u1 = var->c[m][t] + par->alphai[i] * data->s[m][t];

				if (u1 > u0) {
					pa1 = 1 / (1 + exp(u0 - u1));
				}
				else {
					pa1 = 1 - 1 / (1 + exp(u1 - u0));
				}

				
				//				fprintf(outp, "%d\t%d\t%d\t%f\n",m,i,t,pa1);
   
				 /*       if (pa1>0.999999)
				          pa1 = 0.999999;
				        else if (pa1<0.000001)
				          pa1 = 0.000001;
						  */

				if (data->d[m][i][t] == 1) {
					llh->mi[m][i] += log(pa1);
					break;
				}
				else {
					llh->mi[m][i] += log(1 - pa1);
					/*
					if (isnan(llh->mi[m][i])) {
					tid = omp_get_thread_num();
					printf("Running on thread %d\n", tid);
					int waitme;
					printf("\n m = %d, i = %d, llh->mi[m][i] = %f", m, i, llh->mi[m][i]);
					scanf("input var ?%d", &waitme);
					}*/
				}
			}// t
			llh->i[i] += llh->mi[m][i];


			/*
			if (isnan(llh->mi[m][i])) {
				tid = omp_get_thread_num();
				printf("Running on thread %d\n", tid);
			int waitme;
			printf("\n m = %d, i = %d, llh->mi[m][i] = %f", m, i, llh->mi[m][i]);
			scanf("input var ?%d", &waitme);
			}*/

		} // m
	} // i
//	fclose(outp);

	//printf("Outof the pragma settings: Number of threads: %d\n", omp_get_num_threads());

	if (flgi==0)
	{
		llh->all = 0;
		if (flgi==0)
		{
			for (m=1;m<=n_mkt;m++)
			{
				for (i=1;i<=n_person;i++)
				{
					llh->all += llh->mi[m][i];

					/*tid = omp_get_thread_num();
					printf("Running on thread %d\n", tid);
					
					if (isnan(llh->mi[m][i])) {
						//int waitme;
						printf("\n m = %d, i = %d, llh->mi[m][i] = %f", m, i, llh->mi[m][i]);
						//scanf("input var ?%d", &waitme);
					}*/
		 		} // i
			} // m
		}
	}

	free_dmatrix(avf, 1, n_mkt, 1, n_t); 
}



/*******************************************************************************
 Function: kernel_app
 *******************************************************************************/
double kernel_app( const struct PAR *par, const struct DATA *data, const struct PSD *psd, const struct RND *rnd,
	const struct VAR *var, const int cm, const int ci, const int ct, const int flg)
{
  //======================================================================
  // Declare variables
  //======================================================================
  int it=0, k=0;// , i, j;
  double ss=0,w=0,mean_s=0;
  double mx_skern = 0;
  double av = 0;
  //double * tkern = dvector(1, n_pvmh);
  double * skern = dvector(1, n_pvmh);

  //======================================================================
  // Approximate expected value function at (ci,ct,cs) ct is already next period t
  //======================================================================
  if (flg==0)
  { // based on data (for llh construction)
    mean_s = par->rho_p * data->s[cm][ct-1];
  }else
  { // based on random draw (for pseudo-vf computation)
    mean_s = par->rho_p * rnd->s[ct-1][rnd->i];
  }


	for (it=1;it<=psd->itmx;it++)
	{
		if (psd->flgitmx == 0 || it < psd->vmh)
		{
			k = psd->vmh - it;
		}
		else
		{
			k = n_pvmh + psd->vmh - it;
		}

		skern[k] = 0;

		// price-coefficient, thetai[i][1], apply for any cm, but only ci
		skern[k] -= (par->alphai[ci] - psd->alphai[ci][k])
										*(par->alphai[ci] - psd->alphai[ci][k]);

		// gamma[0], common across individuals and markets
		skern[k] -= (par->gamma[0] - psd->gamma[0][k])
										*(par->gamma[0] - psd->gamma[0][k]);
		
		// gamma[1], common across individuals and markets
		skern[k] -= (par->gamma[1] - psd->gamma[1][k])
										*(par->gamma[1] - psd->gamma[1][k]);

		// phi, common across individuals and markets
		skern[k] -= (par->phi - psd->phi[k])*(par->phi - psd->phi[k]);


		// transition density (note random draw common across markets)
		// tkern[k]  = (psd->rs[ct][k] - mean_s)*(psd->rs[ct][k] - mean_s);

		// include transition prob

		//psd->skern[k] = inv_2bwh * psd->skern[k] - 0.5 * psd->tkern[k] / par->var_gamma;
		skern[k] = inv_2bwh * skern[k] - 0.5 * ((psd->rs[ct][k] - mean_s)*(psd->rs[ct][k] - mean_s)) / par->var_p;

		if (it == 1)
			mx_skern = skern[k];
		else
			mx_skern = max(mx_skern, skern[k]);

	} // it loop

	//(*av) = 0;
	for (it=1; it<=psd->itmx; it++)
	{
		if (psd->flgitmx == 0 || it < psd->vmh)
		{
			k = psd->vmh - it;
		}
		else 
		{
			k = n_pvmh + psd->vmh - it;
		}

		w = exp(skern[k] - mx_skern);
		ss += w;
		(av) += w * psd->pvf[cm][ci][ct][k];
		/*
		if (isnan(av)) 
		{
			int waitme;
			printf("\n kernal hit: w = %f, cm = %d, ci=%d, ct=%d, k = %d, psd->pvf[cm][ci][ct][k] = %f, *av = ", w, cm, ci, ct,k, psd->pvf[cm][ci][ct][k], av);
			scanf("input var ?%d", &waitme);
		}*/
	} // it loop

	if (ss > 0)
		(av) /= ss;

	/*
	if (isnan(av)) {
		int tid = omp_get_thread_num();
		printf("Running on thread %d\n", tid);
		int waitme;
		scanf("input var ?%d", &waitme);
	}*/
	//free_dvector(tkern,     1, n_pvmh);
	free_dvector(skern,     1, n_pvmh);

	return av;
}



/*******************************************************************************
 Function: pseudo_vf
 *******************************************************************************/
void pseudo_vf(const struct PAR *par, const struct DATA *data, struct PSD *psd, struct VAR *var, struct RND *rnd)
{
  //======================================================================
  // Declare variables
  //======================================================================
	int m, i, t; // k;
  double av,u0,u1;

  //======================================================================
  // computation of pseudo-value functions by backward induction
	//	for now, do in full, for each m and i
  //======================================================================
	for (m=1;m<=n_mkt;m++)
	{

		// c_value - compute a sequence of var->c[m][i][t=1,...,n_t]
		c_value(par, data, var, m);

		for (i = 1; i <= 1; i++)
		{
			//		for (i=1;i<=n_person;i++){

			//			psd->ivmh = i;

						// t=T
			t = n_t;

			u0 = 0;
			u1 = var->c[m][t] + par->alphai[psd->ivmh] * rnd->s[t][rnd->i];

			if (u1 > u0)
			{
				psd->pvf[m][psd->ivmh][t][psd->vmh] = log(1 + exp(u0 - u1)) + u1 - euler_c;
			}
			else {
				psd->pvf[m][psd->ivmh][t][psd->vmh] = log(1 + exp(u1 - u0)) + u0 - euler_c;
			}


			// t<T
			for (t=n_t-1;t>=2;t--)
			{

				av = kernel_app(par, data, psd, rnd, var, m, psd->ivmh, t+1, 1);

				u0 = par->beta * av;
				u1 = var->c[m][t] + par->alphai[psd->ivmh] * rnd->s[t][rnd->i];
	
				if (u1 > u0) {
					psd->pvf[m][psd->ivmh][t][psd->vmh] = log(1 + exp(u0 - u1)) + u1 - euler_c;
				}
				else {
					psd->pvf[m][psd->ivmh][t][psd->vmh] = log(1 + exp(u1 - u0)) + u0 - euler_c;
				}

			} // t for t<n_t

		} // i

		for (i=1;i<=n_person;i++)
		{
			for (t=2;t<=n_t;t++)
			{
				psd->pvf[m][i][t][psd->vmh] = psd->pvf[m][psd->ivmh][t][psd->vmh];
			}
		}

	} // m

	for (i=1;i<=n_person;i++)
	{
		psd->alphai[i][psd->vmh] = par->alphai[psd->ivmh];
//		psd->alphai[i][psd->vmh] = par->alphai[i];
	}
	psd->gamma[0][psd->vmh] = par->gamma[0];
	psd->gamma[1][psd->vmh] = par->gamma[1];
	psd->phi[psd->vmh]      = par->phi;

	for (t = 2; t <= n_t; t++)
	{
		psd->rs[t][psd->vmh] = rnd->s[t][rnd->i];
	}

	// forward index for random draws 
	if (rnd->i < n_rnd)
	{
		(rnd->i)++;
	}
	else {
		rnd->i = 1;
	}

	// forward index for person (HETERO ESTIMATION)
	if (psd->ivmh < n_person) {
		(psd->ivmh)++;
	}
	else {
		psd->ivmh = 1;
	}

	// forward index for pseudo-values i
	if (psd->vmh == n_pvmh) {
		psd->vmh = 1;
	}
	else {
		(psd->vmh)++;
	}

	if (psd->vmhmx < n_pvmh)
	{
		(psd->vmhmx)++;
	}

}




/*******************************************************************************
 Function: true_vf
 *******************************************************************************/
void true_vf(const struct PAR *par, const struct DATA *data, const struct PSD *psd, struct VAR *var, const struct RND *rnd,
	const int cm, const int ci)
{
	//======================================================================
	// Declare variables
	//======================================================================
	int t, r;
	double av, u0, u1, dtemp1;// dtemp2;

	//======================================================================
	// computation of true value functions by backward induction
	//======================================================================
	//#pragma omp parallel for schedule(dynamic) private(s,d,k,dtemp1,dtemp2,dtemp3,dtemp4,dtemp5)
	// t=T
	t = n_t;
	u0 = 0;
	dtemp1 = var->c[cm][t];
#pragma omp parallel for schedule(dynamic) private(r, u1, u0) // I've checed that this (and the next) parallelization is correct
	for (r = 1; r <= n_rnd; r++)
	{
		u0 = 0;
		u1 = dtemp1 + par->alphai[ci] * rnd->s[t][r];

		if (u1 > u0)
		{
			var->vf[cm][ci][t][r] = log(1 + exp(u0 - u1)) + u1 - euler_c;
		}
		else
		{
			var->vf[cm][ci][t][r] = log(1 + exp(u1 - u0)) + u0 - euler_c;
		}

		//		if (cm==1 && ci==1)
		//			printf("var->vf[1][1][n_t][%d]=%f\n",r,var->vf[1][1][n_t][r]);

		//			if (cm==1)
		//				printf("r=%d, s=%f, vf=%f\n",r,rnd->s[t][r],var->vf[cm][ci][t][r]);
	} // r



	// t<T until t==2
	for (t = n_t - 1; t >= 2; t--)
	{
		dtemp1 = var->c[cm][t];
		// backward induction here
		// WE MUST HAVE var[t+1] to calculate var[t], so we CAN'T put 
		// OpenMP in outer loop!!!
		#pragma omp parallel for schedule(dynamic) private(r, u0, u1, av)
		for (r = 1; r <= n_rnd; r++)
		{

			av = exp_vf(par, data, psd, rnd, var, cm, ci, t + 1, r);

			u0 = par->beta * av;
			u1 = dtemp1 + par->alphai[ci] * rnd->s[t][r];

			if (u1 > u0) {
				var->vf[cm][ci][t][r] = log(1 + exp(u0 - u1)) + u1 - euler_c;
			}
			else {
				var->vf[cm][ci][t][r] = log(1 + exp(u1 - u0)) + u0 - euler_c;
			}

			//			if (t==2 && cm==1)
			//				printf("r=%d, s=%f, vf=%f\n",r,rnd->s[t][r],var->vf[cm][ci][t][r]);

		} // r
	} // t 
}



/*******************************************************************************
 Function: exp_vf
 *******************************************************************************/
double exp_vf(const struct PAR *par, const struct DATA *data, const struct PSD *psd, const struct RND *rnd,
	const struct VAR *var, const int cm, const int ci, const int ct, const int cr)
{
	//======================================================================
	// Declare variables
	//======================================================================
	int r = 0;
	double mean_s = 0, dtemp0 = 0, dtemp3 = 0, dtemp4 = 0;

	//======================================================================
	// Compute expected value function at (ci,ct,cs) ct is already next period t
	//======================================================================
	if (cr == 0) { // based on data (for llh construction)
		mean_s = par->rho_p * data->s[cm][ct - 1];
	}
	else { // based on random draw (for true value function computation)
		mean_s = par->rho_p * rnd->s[ct - 1][cr];
	}

	dtemp0 = dtemp4 = 0;
	for (r = 1; r <= n_rnd; r++)
	{
		dtemp3 = rnd->s[ct][r] - mean_s;
		dtemp3 = exp(-0.5*dtemp3*dtemp3 / par->var_p);
		dtemp4 += dtemp3;
		dtemp0 += var->vf[cm][ci][ct][r] * dtemp3;
	}
	dtemp0 /= dtemp4;

	return dtemp0;
}


/*******************************************************************************

 Function: c_value
 *******************************************************************************/
void c_value(const struct PAR *par, const  struct DATA *data, struct VAR *var,  int cm)
{
  //======================================================================
  // Declare variables
  //======================================================================
	int t;
  
  //======================================================================
  // Compute var.c[cm][ci][t]
  //======================================================================
	// t=1
	var->c[cm][1] = exp(par->gamma[0] + par->gamma[1] * data->q[cm]);
	// t>1
	for (t = 2; t <= n_t; t++)
	{
		var->c[cm][t] = (1 - par->phi) * var->c[cm][t - 1];
	}
}
