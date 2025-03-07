function [solution] = solve_fixed_point(problem,w,para)
%% We use a fixed point approach to solve the discretized problem

% rng(1);
tic

if nargin < 1
    problem = setup;
    problem.convex = true;
end

%% start vector
A = problem.A;
if nargin < 2
	%start vector can be near A as r(A) = Tp(A)/\alpha might be close to 0.
	if (not(exist('w'))) || (size(A,1) ~= size(w,1))
		if problem.convex
			w = A + randn(size(problem.a,1),1)*1e-7;
		else
			w = A + randn(size(problem.a,1),1);
		end
	end
end

if nargin < 3
	para = struct()
end

% Convergence criterion
maxiter = mygetfield(para, 'maxiter', 10000);
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

iter = 1;
[J, dJ, data] = obj( w );
fprintf('%3s %12s %12s %12s (%12s, %2s)\n', 'it', 'obj', 'grad', 'res', 'ls_alpha', 'it');
list_dJ = [];
list_J = [];
list_r = [];

while true
	% The searching direction is just r(w).
	dw = data.res;
	normdw = norm(dw,'inf');
	normdJ = norm(dJ,'inf');

	if iter == 1
		normdw_tol = min( abstol, normdw * reltol );
	end

	fprintf('%3d %e %e %e (%e, %2d)\n', iter, J, normdJ, normdw, ls_alpha, evals);
	list_dJ = [list_dJ; normdJ];
	list_J = [list_J; J];
	list_r = [list_r; normdw];

	if normdw <= normdw_tol || iter > maxiter
		break
	end

	assert( dJ'*dw < 0 );

	%Wolfe strategy to find suitable step size
	[ls_alpha, evals, J, dJ, data] = ls_strong_wolfe(obj, w, dw, ls_alpha, parameters, false, J, dJ);
	w = w + dw * ls_alpha;

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
