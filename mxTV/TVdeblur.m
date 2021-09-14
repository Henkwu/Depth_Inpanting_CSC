function [X,info] = TVdeblur(B,PSF,delta,eps_rel,rho,gamma)
%TVDEBLUR  Total variation image deblurring
% 
% X = TVdeblur(B,PSF,delta)
% [X,info] = TVdeblur(B,PSF,delta,eps_rel,rho,gamma)
%
% This function solves the TV deblurring problem
%
%    min  TV(X)  subject to   || PSF(*)X - B ||_F <= delta
%
% where B is a blurred noisy image, X is the reconstruction, and delta is
% an upper bound for the residual norm.  The TV function is the 1-norm
% of the gradient magnitude, computed via neighbor pixel differences.
% At the image borders, we imposed reflexive boundary conditions for
% the gradient computations.
%
% PSF(*)X is the image X convolved with the doubly symmetric point spread
% function PSF using reflexive boundary conditions.  In the code, the
% blurring matrix that represents PSF is replaced by a rank-deficient
% well-conditioned approximation obtained by neglecting all eigenvalues
% smaller than rho times the largest eigenvalue (default rho = 1e-3).
%
% The parameter delta should be of the same size as the norm of the image
% noise.  If the image is m-times-n, and sigma is the standard deviation
% of the image noise in a pixel, then we recommend to use delta =
% tau*sqrt(m*n)*sigma, where tau is smaller than one, say tau = 0.5.
%
% The parameter gamma is a an upper bound on the norm of the solution's
% component in the subspace corresponding to the neglected eigenvalues.
% The default value is gamma = sqrt(mn)*max(B(:)) which should be
% sufficient for most problems, see info.STATUS='GAMMA_BOUND_TOO_TIGHT'
% below for expert usage.
%
% The function returns an epsilon-optimal solution X, meaning that
% if X* is the exact solution, then our solution X satisfies
%
%    TV(X) - TV(X*) <= epsilon = max(B(:))*m*n*eps_rel,
%
% where eps_rel is a specified relative accuracy (default eps_rel = 1e-2).
% 
% The solution status is returned in the stuct info, with info.STATUS
% having one of the settings
%  'EPSILON-OPTIMAL-SOLUTION': X is an epsilon-optimal solution
%  'NOT-EPSILON-OPTIMAL-SOLUTION': X is not an epsilon-optimal solution
%  'MAXIMUM-NUMBER-OF-ITERATIONS-EXCEEDED': X is not an epsilon-optimal
%     solution when the maximum number of iterations was reached. 
%  'GAMMA_BOUND_TOO_TIGHT': The artificial bound on X corresponding to
%     small eigenvaules was too tight, most likely because the problem was
%     not sufficiently regularized.  Increase tau and/or decrease rho.
%     Alternatively, increase gamma.
%
% See also: TVdenoise, TVinpaint.

% Other fields of info:
%  info.NDEBLUR            Upper bound for the number of iterations.
%  info.ITERATIONS_K       The number of iterations used.	
%  info.EPS_REL_K          The relative accuracy reached.
%  info.CONDITION_MATRIX   The condition number of the rank-reduced matrix. 
%  info.CARD_I             Number of maintained eigenvalues.
%  info.GAMMA              The chosen gamma.
%  info.NORM_COMP_NEG_EIG  The norm bound on the component of neglected
%                          eigenvalues. If GAMMA = NORM_COMP_NEG_EIG then
%                          'GAMMA_BOUND_TOO_TIGHT' is returned in
%                          info.STATUS, see diagnostics above.
%  info.TIME               Time in seconds for the program to run. 

% J. Dahl^1, P.C. Hansen^2, S.H. Jensen^1 & T.L. Jensen^1
% CSI project: (1) Aalborg University, (2) Technical University of Denmark
% April 28, 2009.

if nargin < 3	
    error('Too few input parameters');
elseif size(PSF) ~= size(B)
    error('The PSF must have the same dimensions as B');
elseif nargin == 3
    eps_rel = 1e-2;
    rho = 1e-3;
    gamma = sqrt(numel(B))*norm(B(:),'inf');
elseif nargin == 4
    rho= 1e-3;
    gamma = sqrt(numel(B))*norm(B(:),'inf');
elseif nargin == 5
    gamma = sqrt(numel(B))*norm(B(:),'inf');
end

tic
    
R = max(abs(B(:)));
[m,n] = size(B);
mn = m*n;

% Compute the DCT spectrum of the PSF.
e1 = zeros(m,n); e1(1,1) = 1.0;
A = dctshift(PSF,[ceil(m/2),ceil(n/2)]);
S = dcts2(A)./dcts2(e1);
maxS = max(abs(S(:)));

% Tjeck if delta = 0, in which the solution is given.
if delta==0
    warning('Computing naive solution that may contain Inf and/or NaN')
    X = idcts2(dcts2(B)./S);
    info = info_type_deblur(1,0,0,0,0,0,0,0,toc);
    return;
end

% Normalize the gamma bound.
gamma = gamma/maxS;

% Find indices corresponding to values larger or smaller than max(|S|)*rho.
Ic = int32(find(abs(S(:)) < maxS*rho ));
I  = int32(setdiff(1:(m*n),Ic));

% More parameters for the deblurring algorithm.
epsilon = m*n*R*eps_rel;
mu  = epsilon/mn;
Lmu = 8/mu ;

if isempty(Ic) % The matrix is well conditioned, no rank reduction.
    
    % Compute the implicit bound via the trust region problem
    %    max 0.5 ||xb||_2^2
    %    s.t. || diag(S)*xb - bb||_2 < delta
        
    dct2B = dcts2(B);
    zmax = mxtrp(-S.*dct2B,-S.^2,delta^2);
    xmax = (zmax + dct2B)./S;
    
    D1 = 0.5*norm(xmax)^2;
    N  = int32( ceil(4*sqrt(8*D1*(mn/2))/epsilon) );
    
    % Compute TV solution via C function.
    [X,k,epsilon_k] = tv_deblur(B,S,delta,epsilon,Lmu,mu,N);
    gamma = Inf;
    norm_comp_neg_eig = 0;
    
else % Use rank reduction to avoid the smallest eigenvalues.
    
    % Compute the implicit bound via the trust region problem
    %    max 0.5 ( ||x_I||_2^2 + ||x_Ic||_2^2 )
    %    s.t. || (diag(S)*x - bb)_I||_2<delta
    %         || x_Ic ||_2 < gamma
    
    dct2B = dcts2(B);
    
    zmax = mxtrp(-S(I).*dct2B(I),-S(I).^2,delta^2);
    xmax = (zmax + dct2B(I))./S(I);
    
    D1 = 0.5*norm(xmax)^2 + 0.5*gamma^2;
    N  = int32( ceil(4*sqrt(8*D1*(mn/2))/epsilon) );

    % Compute TV solution via C function.
    [X,k,epsilon_k] = tv_deblur_rr(B,S,I,Ic, ...
                                   delta,gamma,epsilon,Lmu,mu,N);
   
    dct2X = dcts2(X);
    norm_comp_neg_eig = norm(dct2X(Ic));
    
end

if nargout == 2
   
	if k >= N
        info = info_type_deblur(3,N,k,epsilon_k/(R*m*n),max(abs(S(:)))/ ...
                    min(abs(S(:))),length(I),gamma,norm_comp_neg_eig,toc);
	elseif gamma < 1.001 * norm_comp_neg_eig
        info = info_type_deblur(4,N,k,epsilon_k/(R*m*n),max(abs(S(:)))/ ...
                    min(abs(S(:))),length(I),gamma,norm_comp_neg_eig,toc);
    elseif epsilon_k > epsilon
		info = info_type_deblur(2,N,k,epsilon_k/(R*m*n),max(abs(S(:)))/ ...
                    min(abs(S(:))),length(I),gamma,norm_comp_neg_eig,toc);
	else
		info = info_type_deblur(1,N,k,epsilon_k/(R*m*n),max(abs(S(:)))/ ...
                    min(abs(S(:))),length(I),gamma,norm_comp_neg_eig,toc);
	end
end