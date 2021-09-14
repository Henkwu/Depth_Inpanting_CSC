#include <mex.h>
#include "tools.h"
#include "tv_core.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	register double delta, gamma,d,L, mu, eps;
	register double *y,*x,*kf,*epsilon_kf;
	register int *I,*Ic;
	mxArray *Ym,*Im,*Imc;
	mxArray *zp;
	register int maxiter;  
	int m,n,sI,sIc,i;
  
	if(nrhs != 10)
		printf("Should contain 10 input parameters but has %i\n",nrhs); DRAW

	Ym = (mxArray*)prhs[0]; /* Pointer to matrix structure*/
	y = mxGetPr(Ym); /* Pointer to the matrix data*/

	Im = (mxArray*)prhs[1]; 
	I = (int*)mxGetPr(Im);

	Imc = (mxArray*)prhs[2]; 
	Ic = (int*)mxGetPr(Imc); 
	
	zp = (mxArray*)prhs[3];
	delta = (double)(mxGetScalar(zp));

	zp = (mxArray*)prhs[4];
	gamma = (double)(mxGetScalar(zp));

	zp = (mxArray*)prhs[5];
	d = (double)(mxGetScalar(zp));

	zp = (mxArray*)prhs[6];
	eps = (double)(mxGetScalar(zp));

	zp = (mxArray*)prhs[7];
	L = (double)(mxGetScalar(zp));

	zp = (mxArray*)prhs[8];
	mu = (double)(mxGetScalar(zp));

	zp = (mxArray*)prhs[9];
	maxiter = (int)(mxGetScalar(zp));

	m = mxGetM(Ym), n = mxGetN(Ym);
	sI = mxGetM(Im)*mxGetN(Im), sIc = mxGetM(Imc)*mxGetN(Imc);

	/*Allocate memory and assign output pointer*/
	plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL); /*mxReal is our data-type*/
	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);

	/* Get a pointer to the data space in our newly allocated memory */
	x = mxGetPr(plhs[0]);
	kf = mxGetPr(plhs[1]);
	epsilon_kf = mxGetPr(plhs[2]);

	/* Difference between indexing in Matlab and C*/
	for (i=0; i<sIc; i++)
		Ic[i] = Ic[i]-1;

	for (i=0; i<sI; i++)
		I[i] = I[i]-1;

    if(delta ==0)
    	tv_inpaint_core(x,y,I,Ic,gamma,d,eps,L,mu,m,n,sI,sIc,maxiter,kf,epsilon_kf);
    else
        tv_denoise_inpaint_core(x,y,I,Ic,delta,gamma,d,eps,L,mu,m,n,sI,sIc,maxiter,kf,epsilon_kf);
    
  /* Difference between indexing in Matlab and C, revert*/
	for (i=0; i<sIc; i++)
		Ic[i] = Ic[i]+1;

	for (i=0; i<sI; i++)
		I[i] = I[i]+1;
}
