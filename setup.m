function problem = setup(refs, problem)
% sets up all needed data.

%refs       number of refinements mesh will undergo.
%problem    if this is given the existing data (D, F, \Omega, \alpha, ...)
%           is preserved and the full problem structure will be constructed.

if nargin < 2
	problem = struct();
end

%% Define polygon D. D is the domain of u_d
problem = providefield(problem, 'D', 'vertices', ...
            [-1,-1;
              1,-1;
              1, 1;
             -1, 1]...
);

problem.D = convex_polygon_VtoH(problem.D.vertices);

%% Define polygon Omega where y lives.
problem = providefield(problem, 'vertices_Omega', problem.D.vertices);

%% Define polygon F. F is the domain of u. For simplicity F \subset \cl(\Omega)
problem = providefield(problem, 'F', 'vertices', problem.vertices_Omega);
problem.F = convex_polygon_VtoH(problem.F.vertices);

%% Set up the mesh. Points p, edges e and triangle t.
% polygonOmega = @(varargin) polyg(problem.vertices_Omega, varargin{:});
polygonOmega = @(varargin) ourg(problem.vertices_Omega, problem.F.vertices, varargin{:});
[g.p, g.e, g.t] = initmesh(polygonOmega, 'hmax', 2.0);

% Refine the mesh
if nargin < 1
	refs = 4;
end

for ref = 1:refs
	[g.p, g.e, g.t] = refinemesh(polygonOmega, g.p, g.e, g.t);
end

%% Define the discretized F where the a_i live. F is subset of closure of Omega.
a_ind = find(all( problem.F.A * g.p <= problem.F.b )');
a = g.p(:,a_ind)';

% Set up the stiffness and mass matrix.
[K,M] = assema(g.p,g.t,1,1,0);

% number of nodes of Omega.
n = size(g.p,2);

% Extract indices of the interior nodes. We have to be careful in case that
% we added edges of F that are not edges of \Omega to the mesh as these are
% edges inside and would lead to wrong boundary conditions.
outer_mesh_indices = find(any(g.e(6:7,:)==[0;0]));
p_boundary = sort(g.e(1,outer_mesh_indices)');
p_interior = setdiff((1:n),p_boundary);

ni = length(p_interior);

% Submatrices corresponding to inner nodes.
Ki = K(p_interior, p_interior);
Mi = M(p_interior, p_interior);

% Some matrices to convex between vectors corresponding to all nodes, interior nodes or the points a_i.
I = speye(n,n);
I2 = I(:,p_interior);       %interior to all, i.e. I2 : R^ni -> R^n
I3 = I(a_ind,:);            %all to a_ind, i.e. I3 : R^n -> R^aind
I4 = I3*I2;                 %interior to a_ind, i.e. I4 : R^ni -> R^aind

% The desired state
problem = providefield(problem, 'yd_fun', @(x)(cos(pi/2*x(1,:)).*cos(pi/2*x(2,:)).*exp(x(2,:)) / 2));

% The second desired state for the non-convex situation.
problem = providefield(problem, 'yd2_fun', @(x)(0.5*cos(pi/2*x(1,:)+pi).*sin(pi*x(2,:)) / 2));

% Values of desired states at the interior nodes
ydi = problem.yd_fun( g.p(:, p_interior) )';
yd2i = problem.yd2_fun( g.p(:, p_interior) )';

% If problem.convex == false, we have
%   g(u) = .5 || y - yd ||^2,
% otherwise
%   g(u) = .25 || y - yd ||^2 || y - yd2 ||^2.
problem = providefield(problem, 'convex', true);


% Weight parameters
problem = providefield(problem, 'alpha', 1e-3);

%booleans for turning on/off some plots.
problem = providefield(problem,'plot_u',true);
problem = providefield(problem,'plot_residuum',false);
problem = providefield(problem,'plot_p',false);
problem = providefield(problem,'plot_y',false);
problem = providefield(problem,'plot_yd',false);
problem = providefield(problem,'plot_yd2',false);

% Store the data in the problem struct
problem.Ki = Ki;
problem.Mi = Mi;
problem.M = M;
problem.I2 = I2;
problem.I3 = I3;
problem.I4 = I4;
problem.g = g;
problem.ydi = ydi;
problem.yd2i = yd2i;
problem.polygonOmega = polygonOmega;
problem.n = n;
problem.ni = ni;
problem.a = a;
problem.a_ind = a_ind;
problem.p_interior = p_interior;
problem.A = .5*sum(problem.a.^2,2);

end
