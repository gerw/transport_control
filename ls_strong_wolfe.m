function [t, evals, f, g, varargout] = ls_strong_wolfe(fhandle, x, d, alpha, parameters, visualize, f0, g0)
% function alpha = ls_strong_wolfe(fhandle, x, d, alpha, parameters, visualize)
%
% strong Wolfe linesearch with cubic interpolation.
%
%  - x is the starting point
%  - d is the search direction
%  - alpha is the initial step length
%  - parameters contain parameters.sigma and parameters.beta
%  - If the optional parameter visualize is 1, the line seach
%    is visualized.

% Get the parameters
sigma  = parameters.sigma;
tau    = parameters.tau;
gamma   = parameters.gamma;
gamma1 = parameters.gamma1;
gamma2 = parameters.gamma2;

% Keep sure that parameters satisfy the corresponding conditions.
assert( 0 < sigma && sigma < tau && tau < 1 );
assert( 1 < gamma );
assert( 0 < gamma1 && gamma1 <= 1/2 );
assert( 0 < gamma2 && gamma2 <= 1/2 );

% Maximum number of iterations.
maxiter_phase1 = 20;
maxiter_phase2 = 30;

% There is the possibilty to visualize this line-search.
if nargin < 6
	visualize = 0;
end

if visualize
	tried_b = [alpha];
	tried_t = [];
end

% Count the function evaluations.
evals = 0;

if nargin < 8
	[f0,g0] = fhandle(x);
	evals = evals + 1; % Is usually already known, and could be passed to this function.
end

% Keep sure that we have a decent direction.
assert(g0'*d < 0);

% Needed several times.
g0d = g0'*d;
sg0d = sigma * g0d;


% Often, the function values are very noise near the optimum.
% Therefore, all function values which are equal up to f_eps are considered equal
% to increase the numerical stability.
f_eps = (1e-12)*(1+abs(f0));

ft = 0;
gt = zeros(size(g0));
varargout = cell(nargout-4,1);

% Evaluation of auxiliary function psi and its derivative.
function [psi_val,psi_der] = psi(t)
	[ft, gt, varargout{1:end}] = fhandle(x + t*d);
	evals = evals + 1;
	psi_val = ft - (f0 + t*sg0d);
	psi_der = gt'*d - sg0d;
end

% Test the strong wolfe conditions.
function [satisfied] = strong_wolfe_conditions( psi_val, psi_der )
	satisfied = ( psi_val <= f_eps && abs( psi_der + sg0d ) <= tau * abs(g0d) );
end

% Phase 1: Find b satisfying strong Wolfe conditions or psi(b) >= 0 or ( psi(b)<= 0 and psi'(b)>=0 )
k = 1;
a = 0;
psi_a_val = 0;
psi_a_der = (1 - sigma)*g0d;

b = alpha;
[psi_b_val, psi_b_der] = psi(b);
while k < maxiter_phase1 && ...
		not( strong_wolfe_conditions(psi_b_val, psi_b_der) ) && ...
		not( psi_b_val >= f_eps || psi_b_der >= 0 )

	% Increase b
	a = b;
	b = gamma*b;

	% Remember the function values for a.
	psi_a_val = psi_b_val; psi_a_der = psi_b_der;

	% Evaluate function and derivative at b.
	[psi_b_val, psi_b_der] = psi(b);

	if visualize
		tried_b = [tried_b; b];
	end

	k = k + 1;
end

if k == maxiter_phase1
	error('Strong Wolfe line search failed in Phase 1.')
end

k = 0;

if strong_wolfe_conditions( psi_b_val, psi_b_der )
	if visualize
		display('Phase 1 was sufficient.')
	end
	t = b;

	% Fall through the next loop.
	% Do not return here, because of the visualization at the end.
	k = Inf;
else
	if visualize
		if psi_b_val >= f_eps
			display('Phase 1: detected point violating Armijo condition.');
		else
			display('Phase 1: detected point with positive derivative.')
		end
	end
end

% Phase 2: Find point satisfying strong Wolfe conditions.
while k < maxiter_phase2
	% Check that the conditions on the interval are satisfied:
	assert( (psi_a_val <= f_eps) && (psi_a_der < 0) && (psi_b_val >= f_eps || psi_b_der >= 0) )

	if ( psi_b_val > 1e30 )
		% This step seems to be way to large
		t = (a+b)/2;
	elseif (psi_a_val < - f_eps || psi_b_val > f_eps )
		% Function values are not noisy.

		do_cubic = (b - a) > 1e-5;

		if do_cubic
			% [a,b]
			% b-a
			% Use cubic interpolation to compute a good guess.
			A = [1 a a^2 a^3; 0 1 2*a 3*a^2; 1 b b^2 b^3; 0 1 2*b 3*b^2];
			% rcond(A)
			B = [psi_a_val; psi_a_der; psi_b_val; psi_b_der];

			% Solve linear system to obtain the coefficients of the Hermite polynomial.
			X = A\B;

			if abs(X(4)) > 1e-10
				% The polynomial is really cubic.
				if psi_b_der > sigma*abs(g0d)
					% In this case, we can use the minimizer of 
					%      psi_I - sigma*(g0'*d),
					% which corresponds to a minimizer of  f  (and not of  psi)
					X(2) = X(2) + sg0d;
				end

				assert( (4*X(3)^2-12*X(2)*X(4))/(36*X(4)^2) > 0 );

				% The stationary points are:
				t = [-X(3)/(3*X(4)) - sqrt((4*X(3)^2-12*X(2)*X(4))/(36*X(4)^2)),...
					-X(3)/(3*X(4)) + sqrt((4*X(3)^2-12*X(2)*X(4))/(36*X(4)^2))];
			else
				% The polynomial is essentially quadratic.
				% Reinterpolate.
				A = [1 a a^2; 0 1 2*a; 1 b b^2; 0 1 2*b];
				B = [psi_a_val; psi_a_der; psi_b_val; psi_b_der];

				% Solve linear system to obtain the coefficients of the quadratic polynomial.
				X = A\B;

				if psi_b_der > sigma*abs(g0d)
					% In this case, we can use the minimizer of 
					%      psi_I - sigma*(g0'*d),
					% which corresponds to a minimizer of  f  (and not of  psi)
					X(2) = X(2) + sg0d;
				end

				% The stationary point is:
				t = -.5*X(2)/X(3);
			end
		else
			t = (a + b) / 2;
		end
	else
		% In this case the function values are very noisy.
		% Do *not* use them for interpolation,
		% but use a linear interpolation of the derivative.
		t = a - psi_a_der * (b-a) / (psi_b_der - psi_a_der);
	end

	% At least one of our guesses should lie in [a,b].
	% assert(any( a <= t & t <= b ) )
	if not(any( a <= t & t <= b ) )
		% k = maxiter_phase2;
		% break
		t = min(max(t,a),b);
	end

	% The first point has precedence.
	if( a <= t(1) && t(1) <= b )
		t = t(1);
	else
		t = t(2);
	end

	% Clip t to get a guaranteed reduction of the interval.
	t = max(t, a + gamma1*(b-a) );
	t = min(t, b - gamma2*(b-a) );

	% Evaluate psi
	[psi_val, psi_der] = psi(t);

	if visualize
		tried_t = [tried_t; t];
	end

	% Check for the strong Wolfe conditions.
	if strong_wolfe_conditions( psi_val, psi_der )
		break
	else
		% Otherwise, set up new interval.
		if psi_val <= f_eps
			if psi_der < 0
				a = t;
				psi_a_val = psi_val; psi_a_der = psi_der;
			else
				b = t;
				psi_b_val = psi_val; psi_b_der = psi_der;
			end
		else
			b = t;
			psi_b_val = psi_val; psi_b_der = psi_der;
		end
	end

	k = k+1;
end

if visualize
	T = linspace(0,1.1*tried_b(end), 101);
	f = zeros(size(T));
	df = zeros(size(T));
	for i = 1:length(T)
		[f(i),g] = fhandle(x + T(i)*d);
		df(i) = g'*d;
	end

	f_b = zeros(length(tried_b),1);
	df_b = zeros(length(tried_b),1);
	for i = 1:length(tried_b)
		[f_b(i),g] = fhandle(x+tried_b(i)*d);
		df_b(i) = g'*d;
	end

	f_t = zeros(length(tried_t),1);
	df_t = zeros(length(tried_t),1);
	for i = 1:length(tried_t)
		[f_t(i),g] = fhandle(x+tried_t(i)*d);
		df_t(i) = g'*d;
	end

	figure(2); clf; hold on;
	plot( T, f, 'b', T, f0 + sg0d*T, 'g');
	plot(tried_b(1), f_b(1), 'ko', tried_b(2:end), f_b(2:end), 'go');
	plot(tried_t, f_t, 'ro');
	xlabel('t');
	legend('f(x + t d)','f(x) + \sigma t f''(x) d', 'Location', 'Best');

	figure(3); clf; hold on;
	plot( T, df, 'b');
	plot( T([1,end]), tau*abs(g0d)*[1;1], 'g-');
	plot( T([1,end]), -tau*abs(g0d)*[1;1], 'g-');
	plot(tried_b(1), df_b(1), 'ko', tried_b(2:end), df_b(2:end), 'go');
	plot(tried_t, df_t, 'ro');
	xlabel('t');
	legend('f''(x + t d) d', '\pm\tau f''(x) d','Location', 'Best');

	% pause
end

if k == maxiter_phase2
	save('failure.mat', 'x', 'd', 'alpha');
	error('Strong Wolfe line search failed in Phase 2.')
end

f = ft;
g = gt;

end
