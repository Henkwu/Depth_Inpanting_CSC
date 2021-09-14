#include <mex.h>
#include "tools.h"
#include "tv_core.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	register double delta, gamma, L, mu, eps,*kf,*epsilon_kf;
	register double *y,*x,*s;
	register int *J,*Jc;
	mxArray *Ym,*Sm,*Jm,*Jmc;
	mxArray *zp;
	int maxiter,i,m,n,sJ,sJc;  

	if(nrhs != 10)
  		printf("Should contain 10 input parameters but has %i\n",nrhs); DRAW
    
        
	Ym = (mxArray*)prhs[0]; /* Pointer to matrix structure*/
	y = (double*)mxGetPr(Ym); /* Pointer to the matrix data*/

	Sm = (mxArray*)prhs[1];
	s = (double*)mxGetPr(Sm);

	Jm = (mxArray*)prhs[2]; 
	J = (int*)mxGetPr(Jm);

	Jmc = (mxArray*)prhs[3]; 
	Jc = (int*)mxGetPr(Jmc); 

	zp = (mxArray*)prhs[4];
	delta = (double)(mxGetScalar(zp));

	zp = (mxArray*)prhs[5];
	gamma = (double)(mxGetScalar(zp));

	zp = (mxArray*)prhs[6];
	eps = (double)(mxGetScalar(zp));

	zp = (mxArray*)prhs[7];
	L = (double)(mxGetScalar(zp));

	zp = (mxArray*)prhs[8];
	mu = (double)(mxGetScalar(zp));

	zp = (mxArray*)prhs[9];
	maxiter = (int)(mxGetScalar(zp));

    m = mxGetM(Ym), n = mxGetN(Ym);
	sJ = mxGetM(Jm)*mxGetN(Jm), sJc = mxGetM(Jmc)*mxGetN(Jmc);

	/*Allocate memory and assign output pointer*/
	plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL); /*mxReal is our data-type*/
	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);

	/* Get a pointer to the data space in our newly allocated memory */
	x = mxGetPr(plhs[0]);
	kf = mxGetPr(plhs[1]);
	epsilon_kf = mxGetPr(plhs[2]);

  /* Difference between indexing in Matlab and C*/

	for (i=0; i<sJ; i++)
		J[i] = J[i] - 1;

	for (i=0; i<sJc; i++)
	Jc[i] = Jc[i] - 1;

	tv_deblur_rr_core(x,y,s,J,Jc,delta,gamma,eps,L,mu,m,n,sJ,sJc,maxiter,kf,epsilon_kf);

     /* Difference between indexing in C and Matlab*/
	for (i=0; i<sJ; i++)
		J[i] = J[i] + 1;

	for (i=0; i<sJc; i++)
	   Jc[i] = Jc[i] + 1;  
}
