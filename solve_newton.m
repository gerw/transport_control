function [solution] = solve_newton(problem, w, para)
%% We use a Newton method applied to the fixed point equation

% rng(1);
tic

if nargin < 1
	problem = setup;
end

A = problem.A;

%% start vector
if nargin < 2 || size(A,1) ~= size(w,1)
	%start vector can be near A as r(A) = Tp(A)/\alpha might be close to 0.
	if problem.convex
		w = A + randn(size(problem.a,1),1)*1e-7;
	else
		w = randn(size(problem.a,1),1)*10;
	end
end

if nargin < 3
	para = struct()
end

% Convergence criterion
maxiter = mygetfield(para, 'maxiter', 20);
abstol = mygetfield(para, 'abstol', 1e-6);
reltol = mygetfield(para, 'reltol', Inf);


%handle that gives function value and derivative of J at w.
obj = @(w)(eval_objective(problem, w));
% For line search
ls_alpha = 1;
evals = 0;
% Wolfe:
parameters.sigma = 1e-2;
parameters.tau = 0.9;
parameters.gamma = 2.0;
parameters.gamma1 = 0.1;
parameters.gamma2 = 0.1;
% matrices
ni = problem.ni;
K = problem.Ki;
M = problem.Mi;
I4 = problem.I4;



iter = 1;
[J, dJ, data] = obj( w );
fprintf('%3s %12s %12s %12s (%12s, %2s)\n', 'it', 'obj', 'grad', 'res', 'ls_alpha', 'it');
list_dJ = [];
list_J = [];
list_r = [];

while true
	normres = norm(data.res, 'inf');
	normdJ = norm(dJ,'inf');

	if iter == 1
		normres_tol = min( abstol, normres * reltol );
	end

	fprintf('%3d %e %e %e (%e, %2d)\n', iter, J, normdJ, normres, ls_alpha, evals);
	list_dJ = [list_dJ; normdJ];
	list_J = [list_J; J];
	list_r = [list_r; normres];

	if normres <= normres_tol || iter > maxiter
		break
	end



	% if true
	if problem.convex
		% We need direction -r'(w)^{-1} r(w).
		% We linearized -M(y-ydi)+Kp=0, -w+I2/alpha p + A = -r(w) and Ky-u=0 and got
		% S = [-M O1' K; O1 -I I2/problem.alpha; K -Theta(problem.p_interior,:) O2]
		% dywp = S\[zeros(ni,1);-res;zeros(ni,1)]; %[dy, dw, dp]
		% dw = dywp(ni+1:n+ni);
		% It turned out that dw is the Newton direction. However calculation is
		% slow. We modified the system and got the sparse symmetric system

		S = [-M K; K -I4'*data.Theta*I4/problem.alpha];
		assert(issymmetric(S))
		dyp = S\[zeros(ni,1);I4'*(data.Theta*data.res)];
		dp = dyp(ni+1:2*ni);
	else
		% We have to replace the matrix M in S = [-M K; K -I4'*data.Theta*I4/problem.alpha]
		% by the Hessian of g as a function of y. We get
		% (norm1 + norm2) M + M(y-yd2i) * (M(y-ydi))^T + M(y-ydi) * (M(y-yd2i))^T.
		% This matrix is not sparse.
		% We define v1 := Mi*( y - ydi), v2 := Mi*( y - yd2i). Thus
		% S = [-(M(n1+n2) + v1 v2^T + v2 v1^T) K; K -I4'*Theta*I4/problem.alpha].
		% We define d1 := v1+v2, d2 = v1-v2. Thus, v1 v2^T + v2 v1^T = 1/2 * ( d1*d1^T - d2*d2^T ).
		% The first line in S is (-M(n1+n2) - d1 d1^T/2 + d2 d2^T/2)dy Kdp = 0.
		% We define h1 = -d1^T * dy / 2, h2 = d2^T * dy / 2 for a reformulation. Thus for a new S
		% [-M(n1+n2)      K               d1  d2 ] [dy]
		% [    K     -I4'*Theta*I4/alpha  0   0  ] [dp]
		% [    d1^T       0               2   0  ] [h1]
		% [    d2^T       0               0   -2 ] [h2]
		d1 = data.v1 + data.v2;
		d2 = data.v1 - data.v2;
		S = [-M*(data.norm1+data.norm2)                K                    d1               d2; ...
			          K                  -I4'*data.Theta*I4/problem.alpha zeros(ni,1) zeros(ni,1); ...
			          d1'                      zeros(1,ni)                     2                0; ...
			          d2'                      zeros(1,ni)                     0               -2];
		assert(issymmetric(S))
		dyp = S\[zeros(ni,1);I4'*(data.Theta*data.res); 0; 0];
		dp = dyp(ni+1:2*ni);
	end

	%Searching direction
	dw = I4*dp/problem.alpha + data.res;

	if dJ'*dw < -1e-14
		% Sufficient decrease, use Wolfe line search
		[ls_alpha, evals, J, dJ, data] = ls_strong_wolfe(obj, w, dw, 1.0, parameters, false, J, dJ);
		w = w + dw * ls_alpha;
	else
		% No decrease, use a damped step
		[ls_alpha, evals, J, dJ, data] = ls_backtrack(obj, w, dw, 1.0, parameters, J, dJ);
		w = w + dw * ls_alpha;
	end

	iter = iter + 1;
end


%save solution data in structure
solution.time = toc;
solution.list_dJ = list_dJ;
solution.list_J = list_J;
solution.list_r = list_r;
solution.iter = iter;
solution.w = w;
solution.u = data.u;
solution.y = data.y;
solution.p = data.p;

do_plots(problem,solution);
end
