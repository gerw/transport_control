function [x,y] = ourg(vOmega, vF, bs, s)
% This helper function can be used to create
% a domain in the following sense. F ist a subset of \cl(\Omega). We want
% all vertices of F and Omega to be included as well as the segments.
% Example: \Omega = (-1,1)^2, F = [0,1]^2 shall give all 7 points and 8
% segments.

%obtain inequality system for \Omega.
Omega = convex_polygon_VtoH(vOmega);

%number of vertices vF and vOmega.
nFbs = size(vF,1);
nOmegabs = size(vOmega,1);


%% As we have to include all edges of F, we start with these and add the edges of Omega later.

% Number of boundary segments
nbs = size(vF,1);
boundaries = zeros(nbs, 4);
bv = zeros(nbs,2); %vertices corresponding to that boundary edge.
for i = 1:nbs
	j=i+1;
	if i == nbs
		j = 1;
	end
	leng = norm(vF(i,:)-vF(j,:));
	%we only have to decide if the edge of F lies on an edge of \Omega as
	%it will be an edge between inside (1) and outside (0) area, then.
	if any( (Omega.A * vF(i,:)' == Omega.b) & (Omega.A * vF(j,:)' == Omega.b) )
		boundaries(i,:) = [0, leng, 1, 0]; %edge outside.
	else
		boundaries(i,:) = [0, leng, 1, 1]; %edge inside.
	end
	bv(i,:) = [i,j];
end

%% We add the edges of \Omega now.
for k = 1:nOmegabs
	vFonedge = Omega.A(k,:) * vF' == Omega.b(k); %which vF lie on edge k of \Omega?
	sumvFonedge = sum(vFonedge);

	l = k+1;
	if k == nOmegabs
		l = 1;
	end

	if sumvFonedge == 0 %no vF on this edge: we just add edge k.
		nbs = nbs + 1;
		leng = norm(vOmega(k,:) - vOmega(l,:));
		boundaries = [boundaries; [0, leng, 1, 0]];
		bv = [bv; [nFbs+k,nFbs+l]];
	elseif sumvFonedge == 1
		i = find(vFonedge); %vertex vF_i is on edge k.
		if (vF(i,:) == vOmega(k,:)) || (vF(i,:) == vOmega(l,:)) %vertex of F exactly one of the vertices of the edge.
			nbs = nbs + 1;
			leng = norm(vOmega(k,:) - vOmega(l,:));
			boundaries = [boundaries; [0, leng, 1, 0]];
			bv = [bv; [nFbs+k,nFbs+l]];
		else
			%vertex vF_i is in relative interior of edge k.
			%However, only one vF is on edge k. Thus edge k is split in two.
			nbs = nbs+2;
			leng1 = norm(vOmega(k,:) - vF(i,:));
			leng2 = norm(vF(i,:) - vOmega(l,:));
			boundaries = [boundaries; [0,leng1,1,0], [0,leng2,1,0]];
			bv = [bv; [nFbs+k,i]; [i,nFbs+l]];
		end
	elseif sumvFonedge >= 2
		indices = find(vFonedge);
		if length(indices) == 2
			i = indices(1);
			j = indices(2);
			if (i == 1) && (j == nFbs)
				i = nFbs;
				j = 1;
			end
		else
			%should not happen. As F is convex, this means that more than
			%one edge of F is on the boundary of \Omega. This can always
			%be replaced by a single edge.

			if not(vFonedge(1) && vFonedge(end))
				%if not first and last index together on edge k, choose
				%smallest und biggest.
				i = min(find(vFonedge));
				j = max(find(vFonedge));
			else
				%first and last are on edge k. Considering the vector on a
				%circle, all zeros and all ones are next to eachother.

				i = max(find(not(vFonedge))) + 1; %The successor of the last vertex not on edge k.
				j = min(find(not(vFonedge))) - 1; %the predecessor of the first vertex not on edge k.
			end
		end

		if (all(vF(i,:) == vOmega(k,:))) && (all(vF(j,:) == vOmega(l,:))) %vertices vF_i and vF_j coincide with vOmega_k and vOmega_l
			%edge already exists.
		elseif all(vF(i,:) == vOmega(k,:)) %vF_i is vOmega_k but not vF_j vOmega_l.
			%edge from vOmega_k = vF_i to vF_j already exists.
			nbs = nbs+1;
			leng2 = norm(vF(j,:) - vOmega(l,:));
			boundaries = [boundaries; [0,leng2,1,0]];
			bv = [bv; [j,nFbs+l]];
		elseif all(vF(j,:) == vOmega(l,:)) %vF_j is vOmega_l but not vF_i vOmega_k.
			%edge from vF_i to vF_j = vOmega_l already exists.
			nbs = nbs+1;
			leng1 = norm(vOmega(k,:) - vF(i,:));
			boundaries = [boundaries; [0,leng1,1,0]];
			bv = [bv; [nFbs+k,i]];
		else
			%vF_i and vF_j both inside edge k.
			nbs = nbs+2;
			leng1 = norm(vOmega(k,:) - vF(i,:));
			leng3 = norm(vF(j,:) - vOmega(l,:));
			boundaries = [boundaries; [0,leng1,1,0]; [0,leng3,1,0]];
			bv = [bv; [nFbs+k,i]; [j,nFbs+l]];
		end
	end
end



if nargin == 2
	x = nbs;
	return
end

if nargin == 3
	x = boundaries(bs, :)';
	return
end

vertices = [vF;vOmega];

x = zeros(size(bs));
y = zeros(size(bs));
for i = 1:nbs
	ii = find(bs == i);
	if length(ii) > 0
		%vertices for this boundary edge.
		k = bv(i,1);
		l = bv(i,2);

		x(ii) = interp1([boundaries(i,1) boundaries(i,2)], [vertices(k,1), vertices(l,1)], s(ii));
		y(ii) = interp1([boundaries(i,1) boundaries(i,2)], [vertices(k,2), vertices(l,2)], s(ii));
	end
end
end
