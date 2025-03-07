function [x,y] = polyg(vertices, bs, s)
% This helper function can be used to create
% an arbitrary, polygonal domain, cf. lshape.g.

% Number of boundary segments
nbs = size(vertices,1);

% For convenience, double the first vertex
vertices = [vertices; vertices(1,:)];

% Setup the boundary structure, see 'help pdegeom'
boundaries = zeros(nbs, 4);
for i = 1:nbs
	boundaries(i,:) = [0, norm(vertices(i,:)-vertices(i+1,:)), 1, 0];
end

if nargin == 1
	x = nbs;
	return
end

if nargin == 2
	x = boundaries(bs, :)';
	return
end

x = zeros(size(bs));
y = zeros(size(bs));
for i = 1:nbs
	ii = find(bs == i);
	if length(ii) > 0
		x(ii) = interp1([boundaries(i,1) boundaries(i,2)], [vertices(i,1), vertices(i+1,1)], s(ii));
		y(ii) = interp1([boundaries(i,1) boundaries(i,2)], [vertices(i,2), vertices(i+1,2)], s(ii));
	end
end
end
