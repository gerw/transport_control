function S = convex_polygon_VtoH(v)
% function S = convex_polygon_VtoH(v)
%
% Computes H representation from V representation
%
% v          Vertices of polygon given in MATHEMATICALLY POSITIVE ORDER,
%
% S is a struct containing the fields:
% vertices       copy of input
% A,b            polygon can be described by {x | Ax <= b}.
% integralnorm   integral of |x|^2 over D.

n = size(v,1);

S = struct();
S.vertices = v;

v = [v; v(1,:)];

A = zeros(n,2);
b = zeros(n,1);
for i=1:n
	% v2-v1 gives the outer normal vector (v2(2)-v1(2), v1(1)-v2(1)).
	A(i,:) = [v(i+1,2)-v(i,2), v(i,1) - v(i+1,1)];
	b(i) = A(i,:) * v(i,:)';
end




%% calculate integral of |x|^2 over D.
% We use Gauß' integration theorem and parametrizations along the boundary edges.
% \int_D |x|² dx = 1/3 \int_D div(x_1³ x_2³) dx = 1/3 \int_{\partial D} x_1³ n1 + x_2³ n2 ds
% = 1/3 \sum_{edge PQ of D where d=Q-P}
%  \int_0^1 (P_1+td_1)³ n1 + (P_2+td_2)³ n2 \dt
intofnorm = 0;
for i=1:n
	P1 = v(i,1);
	P2 = v(i,2);
	d1 = v(i+1,1)-v(i,1);
	d2 = v(i+1,2)-v(i,2);
	n1 = A(i,1);
	n2 = A(i,2);
	% note that the normal vectors do NOT have length 1, however their
	% length is what we need, i.e. length of edge.
	intofnorm = intofnorm + 1/3 * (n1*(P1^3 + 3*P1^2*d1 /2 + 3*P1*d1^2 /3 + d1^3 /4)+ ...
		n2*(P2^3 + 3*P2^2*d2 /2 + 3*P2*d2^2 /3 + d2^3 /4));
end

S.A = A;
S.b = b;
S.intofnorm = intofnorm;

end
