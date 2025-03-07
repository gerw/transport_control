function [edges, points, Areas] = compute_dominating_regions( problem, w)
% Computes the structure of the sets on which the affine functions attain the supremum

% problem    problem structure,
% problem.a  contains vectors a_i,
% problem.D  contains data regarding polygon D,
% w          weights.

% edges      edges between regions (by index: left region, right region, start point, end point),
% points     coordinates of points where at least three regions touch,
% Areas      areas of the regions,

% D described by its vertices and by the respective system Ax<=b.
vD = problem.D.vertices;
A = problem.D.A;
b = problem.D.b;

% F
a = problem.a;



% We start by computing a convex hull.

X = [a, w];
K = convhulln(X);

% The polyhedron K is given by a triangulation of its facets which are given by
% 3 indices of corners.
% We consider a facet (i,j,k) with vertices bi = (a_i,w_i)^T.
% The OUTER normal vector on this facets is \tilde n = (bk-bi) x (bj-bi).
% We write \tilde n = (normal1, normal2, normal3) and n = (normal1, normal2).
% By orthogonality: \tilde n^T (bk-bi) = 0.
% normal3 \ne 0 follows from fact that ak-ai and aj-ai are not parallel. Thus,
% (-n/normal3)^T (ak-ai) - (wk-wi) = 0.
% Hence, s = -n/normal3 is a point where i,j,k have same value. Still open: is it active?
% K is convex, the projection of bi+\tilde n onto K is bi, hence
% (P(bi+\tilde n) - (bi+\tilde n))^T (x-P(bi+\tilde n)) >= 0 for all x \in K.
% I.e. with bl \in K for all l:
%   -\tilde n^T (bl-bi) >= 0.
%   n^T ai + normal3 wi >= n^T bl + normal3 wl.
% if normal3 is negative we can divide by -normal3.
%   (-n/normal3)^T ai - wi >= (-n/normal3)^T al - wl.
% if normal3 is positive, the sign changes and the point is not active.

points = zeros(size(K,1),2);
I = false(size(K,1),1);
for i = 1:size(K,1)
	P11 = X(K(i,1),:);
	P22 = X(K(i,2),:);
	P33 = X(K(i,3),:);
	% normal = cross(P3 - P1, P2 - P1);
	normal3 = (P33(1) - P11(1))*(P22(2)-P11(2)) - (P33(2) - P11(2))*(P22(1)-P11(1));
	if normal3 < 0
		normal1 = (P33(2) - P11(2))*(P22(3)-P11(3)) - (P33(3) - P11(3))*(P22(2)-P11(2));
		normal2 = (P33(3) - P11(3))*(P22(1)-P11(1)) - (P33(1) - P11(1))*(P22(3)-P11(3));
		points(i,1) = -normal1 / normal3;
		points(i,2) = -normal2 / normal3;
		I(i) = true;
	end
end
% We throw away all edges that have not been marked, i.e. only active points
% remain.
K(~I,:) = [];
points(~I,:) = [];
n_points = size(points,1);

% The points where three regions come together have been computed. We need to
% find the edges between the points.

edges = zeros(n_points*3,4);

% Search for edges
for i = 1:n_points
	idx1 = K(i,1);
	idx2 = K(i,2);
	idx3 = K(i,3);

	if idx1 < idx2
		edges(3*i-2,:) = [idx1, idx2, 0, i];
	else
		edges(3*i-2,:) = [idx2, idx1, i, 0];
	end
	if idx2 < idx3
		edges(3*i-1,:) = [idx2, idx3, 0, i];
	else
		edges(3*i-1,:) = [idx3, idx2, i, 0];
	end
	if idx3 < idx1
		edges(3*i,:) = [idx3, idx1, 0, i];
	else
		edges(3*i,:) = [idx1, idx3, i, 0];
	end
end

% We would like to...
% edges = sortrows(edges);

%% Since edges contains only integers, we can do better:
% edgevector = ((edges(:,1)*n_faces + edges(:,2))*n_points + edges(:,3))*n_points + edges(:,4);
%% There are at most two vectors for each pair of "first two entries"; in one of them, entry 3 is set and entry 4 in the other
n_faces = length(w);
edgevector = (edges(:,1)*n_faces + edges(:,2))*2 + (edges(:,3)>0);
[~,I1] = sort(edgevector);
edges = edges(I1,:);

% Merge edges
i = 1;
j = 2;
while j <= n_points*3
	if edges(i,1) == edges(j,1) && edges(i,2) == edges(j,2)
		edges(i,3) = edges(j,3);

		if j+1 <= n_points*3
			edges(i+1,:) = edges(j+1,:);
			i = i + 1;
		end
		j = j + 2;
	else
		edges(i+1,:) = edges(j,:);
		i = i + 1;
		j = j + 1;
	end
end
n_edges = i;
edges(i+1:end,:) = [];

% Not all of the computed points are inside the region D.
% Psides denotes on which sides of the edges of D the points lie.
% Pi denotes if a point lies inside D.
Psides = A*points' <= b;
Pi = all(Psides);
c = cumsum(Pi);

% Interior points
nPi = sum(Pi);

ndiff = n_points - nPi;


%% Next we have to clip, i,e. outside edges are to be removed. Parts of edges that lie inside D remain.
for i = 1:n_edges
	if edges(i,3) > 0 && edges(i,4) > 0
		% We have a segment
		P1 = points(edges(i,3),:);
		P2 = points(edges(i,4),:);
        
		%is there one of the equations that P1 and P2 violate together?
		commonno = any(not(Psides(:,edges(i,3)) | Psides(:,edges(i,4))));

		P1i = Pi(edges(i,3));
		P2i = Pi(edges(i,4));
        
		% Correct indices such that they will be correct after removal below
		edges(i,3) = c(edges(i,3));
		edges(i,4) = c(edges(i,4));

		if P1i && P2i
			% Both points inside
			continue
		elseif P1i && not(P2i)
			H1 = -A*P1' + b;
			H2 = A*(P2-P1)';
			% We look along the ray P1+t(P2-P1) that crosses the edges of D.
			% We take the minimal (positive) t. Note: H1 >= 0 as P1 inside.
			I = H2 > 0;
			t = min( H1(I) ./ H2(I) );
			n_points = n_points + 1;
			points(n_points,:) = (1-t)*P1 + t*P2;
			K(n_points,1:2) = edges(i,[2 1]);
			edges(i,4) = n_points - ndiff;
		elseif not(P1i) && P2i
			H1 = -A*P2' + b;
			H2 = A*(P1-P2)';
			% We look along the ray P2+t(P1-P2) that crosses the edges of D.
			% We take the minimal (positive) t. Note: H1 >= 0 as P2 inside.
			I = H2 > 0;
			t = min( H1(I) ./ H2(I) );
			n_points = n_points + 1;
			points(n_points,:) = (1-t)*P2 + t*P1;
			K(n_points,1:2) = edges(i,[1 2]);
			edges(i,3) = n_points - ndiff;
		else
			% Both points outside

			% The condition commonno (i.e. P1 and P2 violate one of the
			% conditions of Ax <= b at the same time) is sufficient for not
			% intersecting D. But there are cases that commonno is false and
			% there is no intersection with D anyways. However, numerical
			% examples show having no intersection almost always meant that
			% commonno is true.

			if commonno
				edges(i,3) = 0;
				edges(i,4) = 0;
			else
				% The following part alone is able to decide if the edge
				% crosses D but checking nods beforehand makes it faster.
				H1 = -A*P1' + b;
				H2 = A*(P2-P1)';

				tmin = 0;
				tmax = 1;

				% Need
				% t*H2 <= H1
				% I = H2 > 0;
				% tmax = min(tmax, min( H1(I) ./ H2(I) ));
				% I = H2 < 0;
				% tmin = max(tmin, max( H1(I) ./ H2(I) ));
				% I = H2 == 0;
				% if any(H1(I) < 0)
				% 	tmin = Inf;
				% end
				for j = 1:length(b)
					if H2(j) > 0
						tmax = min(tmax, H1(j)/H2(j));
					elseif H2(j) < 0
						tmin = max(tmin, H1(j)/H2(j));
					elseif H1(j) < 0
						tmin = Inf;
						break
					end
				end

				if tmin < tmax
					n_points = n_points + 1;
					points(n_points,:) = (1-tmin)*P1 + tmin*P2;
					K(n_points,1:2) = edges(i,[1 2]);

					n_points = n_points + 1;
					points(n_points,:) = (1-tmax)*P1 + tmax*P2;
					K(n_points,1:2) = edges(i,[2 1]);

					edges(i,3:4) = [n_points-ndiff-1, n_points-ndiff];
				else
					edges(i,3) = 0;
					edges(i,4) = 0;
				end
			end
		end
	elseif edges(i,3) > 0
		% A ray
		P1 = points(edges(i,3),:);
		F1 = a(edges(i,1),:);
		F2 = a(edges(i,2),:);

		% Correct index such that it will be correct after removal below
		edges(i,3) = c(edges(i,3));

		d = [F1(2) - F2(2), F2(1) - F1(1)];

		H1 = -A*P1' + b;
		H2 = A*d';

		tmin = 0;
		tmax = Inf;

		% Need
		% t*H2 <= H1
		I = H2 > 0;
		tmax = min(tmax, min( H1(I) ./ H2(I) ));
		I = H2 < 0;
		tmin = max(tmin, max( H1(I) ./ H2(I) ));
		I = H2 == 0;
		if any(H1(I) < 0)
			tmin = Inf;
		end
        
		% if the ray pierces through D, we have to add two points. If the
		% start point of the ray is inside, only one is needed.
		if tmin < tmax
			if tmin > 0
				n_points = n_points + 1;
				points(n_points,:) = P1 + tmin*d;
				K(n_points,1:2) = edges(i,[1 2]);
				edges(i,3) = n_points - ndiff;
			end

			n_points = n_points + 1;
			points(n_points,:) = P1 + tmax*d;
			K(n_points,1:2) = edges(i,[2 1]);

			edges(i,4) = n_points - ndiff;
		else
			edges(i,3:4) = 0;
		end
	else
		% A ray
		P2 = points(edges(i,4),:);
		F1 = a(edges(i,1),:);
		F2 = a(edges(i,2),:);

		% Correct index such that it will be correct after removal below
		edges(i,4) = c(edges(i,4));

		d = -[F1(2) - F2(2), F2(1) - F1(1)];

		H1 = -A*P2' + b;
		H2 = A*d';

		tmin = 0;
		tmax = Inf;

		% Need
		% t*H2 <= H1
		I = H2 > 0;
		tmax = min(tmax, min( H1(I) ./ H2(I) ));
		I = H2 < 0;
		tmin = max(tmin, max( H1(I) ./ H2(I) ));
		I = H2 == 0;
		if any(H1(I) < 0)
			tmin = Inf;
		end

		if tmin < tmax
			if tmin > 0
				n_points = n_points + 1;
				points(n_points,:) = P2 + tmin*d;
				K(n_points,1:2) = edges(i,[2 1]);
				edges(i,4) = n_points - ndiff;
			end

			n_points = n_points + 1;
			points(n_points,:) = P2 + tmax*d;
			K(n_points,1:2) = edges(i,[1 2]);

			edges(i,3) = n_points - ndiff;
		else
			edges(i,3:4) = 0;
		end
	end
end

% Remove outside points.
points(~Pi,:) = [];
K(~Pi,:) = [];

% Remove empty edges
edges( all(edges(:,3:4)==0,2),: ) = [];
n_edges = size(edges,1);

%% Add vertices of D
points(end+1:end+size(vD,1),:) = vD;
n_points = size(points,1);

% Boundary points
nPb = n_points - nPi;

% We report the corresponding face
[~,I] = max(a * vD' - w);
K(end+1:end+size(vD,1),1) = I;

%% Add boundary edges
% We calculate the midpoint of D (is inside). We calculate the angles of the
% points on the boundary (vertices of D and new points from clipping). We
% order the angles and get an order of the boundary points, we want to add
% the edges between neighboring boundary points.
midP = sum(vD) / size(vD,1);
phis = atan2( points(nPi+1:end,2)-midP(2), points(nPi+1:end,1)-midP(1) );
[~,I] = sort(phis);

% enlarge array
edges(n_edges + nPb,1) = 0;
for i=1:nPb
	if i < nPb
		ip1 = i + 1;
	else
		ip1 = 1;
	end

	n_edges = n_edges + 1;

	idx1 = K(nPi+I(i), 2);
	if idx1 == 0
		idx1 = K(nPi+I(i), 1);
	end
	idx2 = K(nPi+I(ip1), 1);

	assert( idx1 == idx2 );

	edges(n_edges,:) = [idx1, 0, nPi + I(i), nPi + I(ip1)];
end

if nargout > 2
	% Compute areas
	Areas = zeros(n_faces,1);
	for i = 1:n_edges
		edge = edges(i,:);
		face1 = edge(1);
		face2 = edge(2);
		point1 = edge(3);
		point2 = edge(4);
		P1 = points(point1,:);
		P2 = points(point2,:);

		% Shoelace formula
		A_local = (P1(1)*P2(2) - P2(1)*P1(2))/2;

		if face1
			Areas(face1) = Areas(face1) + A_local;
		end
		if face2
			Areas(face2) = Areas(face2) - A_local;
		end
	end
end
