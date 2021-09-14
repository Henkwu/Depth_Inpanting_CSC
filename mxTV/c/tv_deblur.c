#include <mex.h>
#include "tools.h"
#include "tv_core.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	register double delta, L, mu, eps,*kf,*epsilon_kf;
	register double *y,*x,*s;
	mxArray *Ym,*Sm;
	mxArray *zp;
	register int maxiter;  
	int m,n;
    
	if(nrhs != 7)
		printf("Should contain 7 input parameters but has %i\n",nrhs); DRAW

	Ym = (mxArray*)prhs[0]; /* Pointer to matrix structure*/
	y = mxGetPr(Ym); /* Pointer to the matrix data*/

	Sm = (mxArray*)prhs[1];
	s = mxGetPr(Sm);

	zp = (mxArray*)prhs[2];
	delta = (double)(mxGetScalar(zp));

	zp = (mxArray*)prhs[3];
	eps = (double)(mxGetScalar(zp));

	zp = (mxArray*)prhs[4];
	L = (double)(mxGetScalar(zp));

	zp = (mxArray*)prhs[5];
	mu = (double)(mxGetScalar(zp));

	zp = (mxArray*)prhs[6];
	maxiter = (int)(mxGetScalar(zp));

	m = mxGetM(Ym), n = mxGetN(Ym);

	/*Allocate memory and assign output pointer*/
	plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL); /*mxReal is our data-type*/
	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);

	/* Get a pointer to the data space in our newly allocated memory */
	x = mxGetPr(plhs[0]);
	kf = mxGetPr(plhs[1]);
	epsilon_kf = mxGetPr(plhs[2]);

	tv_deblur_core(x,y,s,delta,eps,L,mu,m,n,maxiter,kf,epsilon_kf);
}
