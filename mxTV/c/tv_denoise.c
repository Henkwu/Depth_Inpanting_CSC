#include <mex.h>
#include "tools.h"
#include "tv_core.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	register double delta, L, mu, eps;
	register double *y,*x,*epsilon_kf,*kf;
	mxArray *Ym;
	mxArray *zp;
	register int maxiter;  
	int type,m,n;

	if(nrhs != 7)
		printf("Should contain 7 input parameters but has %i\n",nrhs); DRAW

	Ym = (mxArray*) prhs[0]; /* Pointer to matrix structure*/
	y = mxGetPr(Ym); /* Pointer to the matrix data*/
	
	zp = (mxArray*) prhs[1];
	delta = (double)(mxGetScalar(zp));

	zp = (mxArray*) prhs[2];
	eps = (double)(mxGetScalar(zp));

	zp = (mxArray*) prhs[3];
	L = (double)(mxGetScalar(zp));

	zp = (mxArray*) prhs[4];
	mu = (double)(mxGetScalar(zp));

	zp = (mxArray*) prhs[5];
	maxiter = (int)(mxGetScalar(zp));

	zp = (mxArray*) prhs[6];
	type = (int)(mxGetScalar(zp));

	m = mxGetM(Ym), n = mxGetN(Ym);

	/*Allocate memory and assign output pointer*/
	plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL); /*mxReal is our data-type*/
	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);

	/* Get a pointer to the data space in our newly allocated memory */
	x = mxGetPr(plhs[0]);
	kf = mxGetPr(plhs[1]);
	epsilon_kf = mxGetPr(plhs[2]);

	if(type==1){
		tv_denoise_core(x,y,delta,eps,L,mu,m,n,maxiter,kf,epsilon_kf);
	}
	else{
		tv_denoise_core_org(x,y,delta,eps,L,mu,m,n,maxiter,kf,epsilon_kf);
	}

}
