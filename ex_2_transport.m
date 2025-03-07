% Helper script to plot solution for some discretization and also the
% transport realized by \nabla \phi_w.
% NOTE: If size of window of figure containing the arrows is changed, arrows
% are moved.

clear;

%% create problem
refs = 2;
%define basic setting of problems.
basic_problem.D.vertices = ...
            [-1,-1;
              1,-1;
              1, 1;
             -1, 1];
basic_problem.vertices_Omega = basic_problem.D.vertices;
basic_problem.F.vertices = basic_problem.D.vertices;

problem = setup(refs,basic_problem);
problem.alpha = 1e-3;
solution = solve_newton(problem);



%% extract data.
w = solution.w;
a = problem.a;
u = solution.u;
g = problem.g;

%% plot solution.
fig1 = figure(1); clf
pdeplot(g.p, g.e, g.t, 'xydata', u, 'zdata', u, 'mesh', 'on');
% title('control $\bar u$', 'interpreter', 'latex', 'FontSize', 28)
set(gca, 'LooseInset', get(gca, 'TightInset'));
set(gca,'XTick',[], 'YTick', [],'ZTick', [])
colormap(parula);
colorbar('off');
xlabel('$x_1$','interpreter', 'latex','FontSize',80,'Position', [mean(xlim) min(ylim) min(u)]);
ylabel('$x_2$','interpreter', 'latex','FontSize',80,'Position', [min(xlim) mean(ylim) min(u)]);

%get and plot edges and points of regions for optimal w.
[edges,points,Areas] = compute_dominating_regions( problem, w );
fig2 = figure(2);
clf; hold on; axis equal;
for i=1:size(edges,1)
	plot( points(edges(i,3:4),1), points(edges(i,3:4),2), '-b',LineWidth=2);
end

%omit a_i where region has no size.
a = (a ./ Areas) .* Areas;

%% compute centers of mass as representation points of the regions.
n_areas = size(Areas,1);
rep = zeros(n_areas,2);

n_edges = size(edges,1);
for i = 1:n_edges
	left = edges(i,1);
	right = edges(i,2);
	P1 = points(edges(i,3),:);
	P2 = points(edges(i,4),:);
	local_dA = P1(1)*P2(2) - P2(1)*P1(2);
	if left > 0 %Not outer region.
		rep(left,:) = rep(left,:) + (P1+P2)*local_dA;
	end
	if right > 0 %Not outer region.
		rep(right,:) = rep(right,:) - (P1+P2)*local_dA;
	end
end
%every point has been counted twice and we divide to get the "middle point".
rep = rep ./ (6*Areas);

%plotting red crosses to mark the a_i
plot(a(:,1), a(:,2), 'rx', 'Markersize',16,'LineWidth',2);
%black circles for representatives.
% plot(rep(:,1), rep(:,2), 'ko');

%plot arrows.
xapf = @(x,pos,xl) pos(3)*(x-min(xl))/diff(xl)+pos(1);                % 'x' Annotation Position Function
yapf = @(y,pos,yl) pos(4)*(y-min(yl))/diff(yl)+pos(2);                % 'y' Annotation Position Function
xl = xlim;
yl = ylim;
pos = gca().Position;
for i=1:n_areas
	if Areas(i) > 0
		annotation('arrow', xapf([rep(i,1) a(i,1)],pos,xl), yapf([rep(i,2) a(i,2)],pos,yl),'LineWidth',2)
	end
end

%change shown area
xlim([-1 1]);
ylim([-1 1]);
xlabel('$x_1$','interpreter', 'latex','FontSize',56,'Position', [mean(xlim) min(ylim) min(u)]);
ylabel('$x_2$','interpreter', 'latex','FontSize',56,'Position', [min(xlim) mean(ylim) min(u)]);
set(gca,'XTick',[], 'YTick', [],'ZTick', [])

% title('Transport of the regions'' masses to the $a_i$','interpreter','latex','FontSize',28);



%save files.
exportgraphics(fig1,"Transport_u_for_" + refs + "_refinements.pdf");
exportgraphics(fig2,"Transport_map_for_" + refs + "_refinements.pdf")
