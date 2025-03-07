function [J, dJ, data] = eval_objective( problem, w )
%% Calculates function value J, derivative dJ and other needed data.

% problem    structure that contains needed data like paramters and matrices,
% w          weights.

% J             function value of target functional,
% dJ            its derivative,
% data.edges    edges between regions,
% data.points   points where at least three regions touch,
% data.u        control values, i.e. the areas of the regions,
% data.Theta    derivative of u w.r.t. w,
% data.y        solution of Ki y = u        (on inner nodes)
% data.p        solution of Ki p = Mi(y-yd) (on inner nodes).
% data.W        Wasserstein term.
% data.res      residual r(w) = A-w+I4*p



% Extract variables
Ki = problem.Ki;                    % stiffness matrix on inner nodes,
Mi = problem.Mi;                    % mass matrix on inner nodes,
ydi = problem.ydi;                  % discretized y_d on inner nodes
a = problem.a;                      % n x 2. Contains vectors a_i.
alpha = problem.alpha;              % parameter in front of Wasserstein term in target functional.
I4 = problem.I4;                    % interior to a_ind, i.e. I4 : R^ni -> R^aind
A = problem.A;                      % contains entries 1/2*|a_i|^2.

[edges, points] = compute_dominating_regions( problem, w );

n_faces = length(w);
n_edges = size(edges,1);

% We have to integrate |x - \nabla u(x)|^2 = |x|^2 + |\nabla u(x)|^2 - 2\nabla u(x)*x over D
% The first addend is computed in 'setup'
W = problem.D.intofnorm;

% Compute areas and their derivatives.
Areas = zeros(n_faces,1);
% Theta will have the diagonal entries (n_faces) and two more for each edge.
I = zeros(2*n_edges + n_faces,1);
J = zeros(2*n_edges + n_faces,1);
S = zeros(2*n_edges + n_faces,1);
diagonal = zeros(n_faces,1);
idx = n_faces;
for i = 1:n_edges
	face1 = edges(i,1);
	face2 = edges(i,2);
	point1 = edges(i,3);
	point2 = edges(i,4);

	% Coordinates of the end points of the edge between dominating regions.
	p11 = points(point1,1);
	p12 = points(point1,2);
	p21 = points(point2,1);
	p22 = points(point2,2);

	% sum of of vertices of edge
	ps1 = p11 + p21;
	ps2 = p12 + p22;

	% Shoelace formula
	Area_local = (p11*p22 - p21*p12)/2;

	% We add the areas via the edges
	if face1
		% The gradient on face1
		a_l_1 = a(face1,1);
		a_l_2 = a(face1,2);
		a_l_normQ = (a_l_1^2+a_l_2^2);

		Areas(face1) = Areas(face1) + Area_local;

		% Contribution of ||\nabla u(x)||^2 = ||a_l||^2
		W = W + a_l_normQ*Area_local;
	end
	if face2
		% The gradient on face2
		a_r_1 = a(face2,1);
		a_r_2 = a(face2,2);
		a_r_normQ = (a_r_1^2+a_r_2^2);

		Areas(face2) = Areas(face2) - Area_local;

		% Contribution of ||\nabla u(x)||^2
		W = W - a_r_normQ*Area_local;
	end

	% for inner edges we consider the derivative
	if face1 && face2
		% Coordinates of the edge vector
		dp1 = p11 - p21;
		dp2 = p12 - p22;
		% Squared length of edge
		leng = dp1^2 + dp2^2;
		% Squared norm of a_l - a_r
		a_lr_normQ = (a_l_1 - a_r_1)^2 + (a_l_2 - a_r_2)^2;
		% Used for the derivative of area
		dArea_local = sqrt(leng/a_lr_normQ);

		% Prepare matrix Theta
		% Non-diagonal-entries (i,j) are just edge_length / norm(a_i - a_j).
		I(idx+1) = face1;
		I(idx+2) = face2;
		J(idx+1) = face2;
		J(idx+2) = face1;
		S(idx+1) = dArea_local;
		S(idx+2) = dArea_local;
		idx = idx + 2;

		% Diagonal entries are minus sum of same terms for all neighbors.
		diagonal(face1) = diagonal(face1) - dArea_local;
		diagonal(face2) = diagonal(face2) - dArea_local;

		% Contribution of -2\nabla u(x)*x
		W = W - 4/3*(w(face1) - w(face2))*Area_local;
	elseif face1 % face 1 and not face2
		% Contribution of -2\nabla u(x)*x
		W = W - 2/3*(a_l_1*ps1+a_l_2*ps2)*Area_local;
	else % (not face1) and face 2
		% Contribution of -2\nabla u(x)*x
		W = W + 2/3*(a_r_1*ps1+a_r_2*ps2)*Area_local;
	end
end

% We reserved first n_faces entries for diagonal elements.
I(1:n_faces) = 1:n_faces;
J(1:n_faces) = 1:n_faces;
S(1:n_faces) = diagonal;
Theta = sparse(I(1:idx),J(1:idx),S(1:idx),n_faces,n_faces);

% We prolong the vector from the controlled a_ind to all interior nodes of Omega.
u = problem.I3'*Areas;

% Solve state equation.
y = Ki\( problem.I2'*u );


if problem.convex
	% Eval objective
	% J = 1/2 * |y - yd|^2 + alpha/2 W_2.
	J = .5*((y - ydi)'*(Mi*(y - ydi))) + alpha / 2 * W;
	% Solve adjoint equation
	p = Ki\(Mi*( y - ydi));
else
	yd2i = problem.yd2i;
	% j(u) = 1/4 * |y-ydi|^2 * |y-yd2i|^2, Ky = u.
	data.v1 = Mi* ( y - ydi);
	data.v2 = Mi* ( y - yd2i);
	data.norm1 = .5*(y - ydi)'*data.v1;
	data.norm2 = .5*(y - yd2i)'*data.v2;
    
	% Eval objective
	% J = 1/2 * |y - yd|^2 * 1/2 * |y - yd2|^2 + alpha/2 W_2.
	J = data.norm1*data.norm2 + alpha / 2 * W;
	% Solve adjoint equation
	p = Ki\(   ( data.norm1 * data.v2 + data.norm2 * data.v1 )      );
end

%residual and derivative of J.
res = A - w + I4*p/alpha;
dJ = alpha * Theta * res;


%put calculated data at point w into data structure.
data.edges = edges;
data.points = points;
data.u = u;
data.Theta = Theta;
data.y = y;
data.p = p;
data.W = W;
data.res = res;
