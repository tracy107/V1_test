/*******************************************************************************

  This program simulates the data from a very simple dynamic model 
	with finite time horizon.

	Consumers either choose to buy or don't buy

	if consumers buy

	First Version: Sep 19, 2012
  Updated: Jan 22, 2013
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <float.h>
#include <time.h>
#include <sys/time.h>
#include "tools.h"
#include "sim.h"
#include <omp.h>


int
  n_param    =     8,  // # all parameters in param_input_sim.txt
  n_mkt      =    30,  // # markets
	o_n_t      =    30,
  n_t        =   100,  // total # periods (terminal)
  n_person   =   100,  // # simulated consumers
	n_theta    =     1,  // # parameters to be estimated
  n_gamma    =     2,  // # param for state transition density
	n_rnd      =   100;  // total # of grid points PER PRODUCT, be careful of setting if it n_product is large as it's gonna be n_rnd^n_product


const int
	euler_c    = 0.5772156649;

int main(void)
{
  //======================================================================
  // Declare variables
  //======================================================================
	int i, j, t, s, m, r, k, itemp1;
	long seed=-123456;
	double dtemp0, dtemp1, dtemp2, dtemp3, dtemp4, u0, u1, maxs2, mins2, llh;

	double *param;

  struct DATA data;     // store data
  struct PAR par;       // store parameters
	struct VAR var;       // store variables
	struct RND rnd;       // store random numbers

  char param_name[n_param+1][50];
  
  FILE *inp, *outp;

  //======================================================================
  // Memory allocation
  //======================================================================
  param       = dvector(1, n_param);
	data.s      = dmatrix(1, n_mkt, 1, n_t);
	data.q      = dvector(1, n_mkt);
	par.gamma   = dvector(0, n_gamma);
  par.alphai  = dvector(1, n_person);
	rnd.nu      = dmatrix(1, n_person, 1, n_theta);
	var.s       = dmatrix(1, n_t, 1, n_rnd);
  var.c       = dmatrix(1, n_mkt, 1, n_t);
	var.cp      = darray3d(1, n_mkt, 1, n_person, 1, n_t);
	var.vf      = darray4d(1, n_mkt, 1, n_person, 2, n_t, 1, n_rnd);

  //============================================================================
  // Read in the starting parameter values
  //============================================================================
  printf("before reading starting parameter values.\n");

  inp = fopen("param_input_sim.txt","r");
  i=1;

  while (i<=n_param){
    k = fscanf(inp, "%lf%s", &param[i], param_name[i]);
    printf("param[%d] = %f   %s.\n", i, param[i], param_name[i]);
    i++;
  }
	fclose(inp);

  par.beta      = param[1];
  par.gamma[0]  = param[2];
  par.gamma[1]  = param[3];
	par.alpha     = param[4];
  par.sig_alpha = param[5];
  par.phi       = param[6];
  par.rho_p     = param[7];
  par.sd_p      = param[8];
  
  printf("after reading starting parameter values.\n");

  par.var_p = par.sd_p * par.sd_p;

  //============================================================================
  // Data generation [observed states]
  //============================================================================
  maxs2=-DBL_MAX;
  mins2= DBL_MAX;
  for (m=1;m<=n_mkt;m++){
    for (t=1;t<=n_t;t++){
      if (t==1){
        data.s[m][1] = 8 + par.sd_p * gasdev(&seed);
      }else{
        data.s[m][t] = par.rho_p * data.s[m][t-1] + par.sd_p * gasdev(&seed);
      }
      maxs2 = DMAX(maxs2, data.s[m][t]);
      mins2 = DMIN(mins2, data.s[m][t]);
      //printf("data.s[%d][%d]=%f\n",m,t,data.s[m][t]);
    } // t
  } // m

	printf("after generating observed states\n");

  //============================================================================
  // Random variable generation
  //============================================================================
	// consumers
  outp = fopen("nu.txt","w");
	for (i=1;i<=n_person;i++){
	    fprintf(outp, "%d", i);
    for (j=1;j<=n_theta;j++){
  		rnd.nu[i][j] = gasdev(&seed);
	    fprintf(outp, "\t%f", rnd.nu[i][j]);
      //printf("rnd.nu[%d][%d]=%f\n",i,j,rnd.nu[i][j]);
    }
	  fprintf(outp, "\n");
	}
  fclose(outp);

  // alpha
  for (i=1;i<=n_person;i++){
    par.alphai[i] = par.alpha + par.sig_alpha * rnd.nu[i][1];
    //printf("alphai[%d]=%f",i,par.alphai[i]);
  }

  // c
  for (m=1;m<=n_mkt;m++){
		data.q[m] = gasdev(&seed);
    var.c[m][1] = exp(par.gamma[0] + par.gamma[1]*data.q[m]);
    for (t=2;t<=n_t;t++){
      var.c[m][t] = (1-par.phi) * var.c[m][t-1];
    }
  }

  //============================================================================
  // Value function conputation
	//	integration is done not by deterministic grid but by random grid
	//============================================================================
	// initialization
	outp = fopen("state.txt","w+");
  for (r=1;r<=n_rnd;r++){
    for (t=1;t<=n_t;t++){
			fprintf(outp, "%d\t%d",r,t);
			// same states for all t (this may create a problem in pseudo-estimation
//			if (t==1)
				var.s[t][r] = (maxs2-mins2) * ran1(&seed) + mins2;
//			else
//				var.s[t][r] = var.s[t-1][r];

			fprintf(outp, "\t%.10f", var.s[t][r]);
			fprintf(outp, "\n");
		} // r
	} // t
	fclose(outp);
  
	// value function computation by backward induction


	outp = fopen("tvf.txt", "w");
  //#pragma omp parallel for schedule(dynamic) private(s,d,k,dtemp1,dtemp2,dtemp3,dtemp4,dtemp5)
  for (m=1;m<=n_mkt;m++){
		printf("m=%d\n",m);
    for (i=1;i<=n_person;i++){
      // t=T
      t=n_t;
      //printf("m=%d, i=%d, t=%d\n",m,i,t);
      //#pragma omp parallel for schedule(dynamic) private(dtemp0,dtemp1,k) // I've checked that this (and the next) parallelization is correct
      for (r=1;r<=n_rnd;r++){ // grid for consumption

        dtemp0 = 0;
        dtemp1 = var.c[m][t] + par.alphai[i] * var.s[t][r];

        if (dtemp1>dtemp0)
          var.vf[m][i][t][r] = log(1+exp(dtemp0-dtemp1))+dtemp1 - euler_c;
        else
          var.vf[m][i][t][r] = log(1+exp(dtemp1-dtemp0))+dtemp0 - euler_c;

				if (i==1)
					fprintf(outp, "%d\t%d\t%d\t%f\t%f\n",m,t,r,var.s[t][r],var.vf[m][i][t][r]);

      } // r
    
      // t<T until t==2
      for (t=n_t-1;t>=2;t--){
        //printf("m=%d, i=%d, t=%d\n",m,i,t);
        //#pragma omp parallel for schedule(dynamic) private(dtemp0,dtemp1,dtemp3,dtemp4,k)
        for (r=1;r<=n_rnd;r++){

          dtemp0=dtemp1=dtemp4=0;
          for (k=1;k<=n_rnd;k++){
            dtemp3 = var.s[t+1][k] - par.rho_p*var.s[t][r];
            dtemp3 = exp(-0.5*dtemp3*dtemp3/par.var_p);
            dtemp4 += dtemp3;
            dtemp0 += var.vf[m][i][t+1][k] * dtemp3;
          }
          dtemp0 /= dtemp4;

          dtemp0 *= par.beta;
          dtemp1 = var.c[m][t] + par.alphai[i] * var.s[t][r];

          if (dtemp1>dtemp0)
            var.vf[m][i][t][r] = log(1+exp(dtemp0-dtemp1))+dtemp1 - euler_c;
          else
            var.vf[m][i][t][r] = log(1+exp(dtemp1-dtemp0))+dtemp0 - euler_c;

					if (i==1)
						fprintf(outp, "%d\t%d\t%d\t%f\t%f\n",m,t,r,var.s[t][r],var.vf[m][i][t][r]);


        } // r
      } // t
     
    } // i

  } // m
	fclose(outp);

  printf("after value funciton computation.\n");

  //============================================================================
  // Data generation
	//============================================================================
	outp = fopen("data.txt","w+");
	llh = 0;
  
  for (m=1;m<=n_mkt;m++){

    for (i=1;i<=n_person;i++){
  
      for (t=1;t<=o_n_t;t++){

        if (t==n_t){
          dtemp0 = 0;
        }else{
          dtemp0=dtemp1=dtemp4=0;
          for (r=1;r<=n_rnd;r++){
            dtemp3 = var.s[t+1][r] - par.rho_p*data.s[m][t];
            dtemp3 = exp(-0.5*dtemp3*dtemp3/par.var_p);
            dtemp4 += dtemp3;
            dtemp0 += var.vf[m][i][t+1][r] * dtemp3;
          }
          dtemp0 /= dtemp4;
        }

        u0 = par.beta * dtemp0;
        u1 = var.c[m][t] + par.alphai[i] * data.s[m][t];
        
			if (u1>u0)
				dtemp3 = 1/(1+exp(u0-u1));
			else
				dtemp3 = 1 - 1/(1+exp(u1-u0));

			if (ran1(&seed) < dtemp3){
				itemp1 = 1;
				llh += log(dtemp3);
			}else{
				itemp1 = 0;
				llh += log(1-dtemp3);
			}

      fprintf(outp, "%d\t%d\t%d\t%d\t%.10f\t%.10f\t%f\t%f\t%f\t%f\t%f\t%f\n", 
											m, i, t, itemp1, data.s[m][t], data.q[m], dtemp0, dtemp3, par.alphai[i],var.c[m][t],u0,u1);

			if (itemp1==1)
				break;

      } // t

    } // i

  } // m

	printf("llh = %f\n",llh);

	// simply compute choce prob for all periods but for person 1
	// the above can't be used because if all consumers adopt, we don't compute choice prob
	// for the remaining periods
	i = 1;
  for (m=1;m<=n_mkt;m++){
    for (t=1;t<=n_t;t++){

      if (t==n_t){
        dtemp0 = 0;
      }else{
        dtemp0=dtemp1=dtemp4=0;
        for (r=1;r<=n_rnd;r++){
          dtemp3 = var.s[t+1][r] - par.rho_p*data.s[m][t];
          dtemp3 = exp(-0.5*dtemp3*dtemp3/par.var_p);
          dtemp4 += dtemp3;
          dtemp0 += var.vf[m][i][t+1][r] * dtemp3;
        }
        dtemp0 /= dtemp4;
      }

      u0 = par.beta * dtemp0;
      u1 = var.c[m][t] + par.alphai[i] * data.s[m][t];
      
			if (u1>u0)
				dtemp3 = 1/(1+exp(u0-u1));
			else
				dtemp3 = 1 - 1/(1+exp(u1-u0));

			var.cp[m][i][t] = dtemp3;

    } // t
  } // m


	outp = fopen("data_check.txt","w+");
	for (t=1;t<=n_t;t++){
		fprintf(outp, "%d",t);
		for (m=1;m<=n_mkt;m++)
			fprintf(outp, "\t%f\t%f\t%f",data.s[m][t],var.c[m][t],var.cp[m][1][t]);
		fprintf(outp, "\n");
	}
	fclose(outp);


	printf("data generation done!\n");

  //======================================================================
  // Free up Memory
  //======================================================================
  free_dvector(param,      1, n_param);
	free_dmatrix(data.s,     1, n_mkt, 1, n_t);
	free_dvector(data.q,     1, n_mkt);
	free_dvector(par.gamma,  0, n_gamma);
  free_dvector(par.alphai, 1, n_person);
	free_dmatrix(rnd.nu,     1, n_person, 1, n_theta);
	free_dmatrix(var.s,      1, n_t,    1, n_rnd);
  free_dmatrix(var.c,      1, n_mkt, 1, n_t);
  free_darray3d(var.cp,    1, n_mkt, 1, n_person, 1, n_t);
	free_darray4d(var.vf,    1, n_mkt, 1, n_person, 2, n_t, 1, n_rnd);
    
  
	return 0;
}
