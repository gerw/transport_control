function [t, evals, f, g, varargout] = ls_backtrack(fhandle, x, d, alpha, parameters, f0, g0)
% function [t, evals, f, g, varargout] = ls_backtrack(fhandle, x, d, alpha, parameters, f0, g0)
%
% simple backtracking line search
%
%  - x is the starting point
%  - d is the search direction
%  - alpha is the initial step length
%  - parameters contain parameters.sigma and parameters.beta

% Get the parameter
gamma   = parameters.gamma;

% Keep sure that parameter satisfy the corresponding conditions.
assert( 1 < gamma );

% Maximum number of iterations.
maxiter = 20;

% Count the function evaluations.
evals = 0;

if nargin < 7
	[f0,g0] = fhandle(x);
	evals = evals + 1; % Is usually already known, and could be passed to this function.
end

% Keep sure that we have a weakly decent direction.
assert(g0'*d <= 1e-14);

varargout = cell(nargout-4,1);

t = alpha;
[ft, gt, varargout{1:end}] = fhandle(x + t*d);
evals = evals + 1;

k = 0;
while k < maxiter && not( ft <= f0 + 1e-8 )
	t = t / gamma;

	[ft, gt, varargout{1:end}] = fhandle(x + t*d);
	evals = evals + 1;

	k = k + 1;
end

if k == maxiter
	error('Backtracking line search failed.')
end

f = ft;
g = gt;

end
