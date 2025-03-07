%%helper script to solve problem for some number of refinements, plot it
%%and save it afterwards. We use solve_newton as it is the fastest.
clear;


refs = 5;
problem.D.vertices = ...
            [-1,-1;
              1,-1;
              1, 1;
             -1, 1];
problem.vertices_Omega = problem.D.vertices;
problem.F.vertices = problem.D.vertices;

problem = setup(refs,problem);
problem.alpha = 1e-3;

solution = solve_newton(problem);


%% plot
u = solution.u;
g = problem.g;
fig1 = figure(1); clf
pdeplot(g.p, g.e, g.t, 'xydata', u, 'zdata', u, 'mesh', 'on');
% title('control $\bar u$', 'interpreter', 'latex','FontSize',32)
xlabel('$x_1$','interpreter', 'latex','FontSize',80,'Position', [mean(xlim) min(ylim) min(u)]);
ylabel('$x_2$','interpreter', 'latex','FontSize',80,'Position', [min(xlim) mean(ylim) min(u)]);
set(gca,'XTick',[], 'YTick', [],'ZTick', [])
set(gca, 'LooseInset', get(gca, 'TightInset'));
colormap parula
colorbar('off');


Y = problem.I2*(solution.y - problem.ydi);
fig2 = figure(2); clf
pdeplot(g.p, g.e, g.t, 'xydata', Y, 'zdata', Y, 'mesh', 'on');
% title('residual $y - y_{d,i}$', 'interpreter', 'latex','FontSize',32)
xlabel('$x_1$','interpreter', 'latex','FontSize',80,'Position', [mean(xlim) min(ylim) min(Y)]);
ylabel('$x_2$','interpreter', 'latex','FontSize',80,'Position', [min(xlim) mean(ylim) min(Y)]);
set(gca,'XTick',[], 'YTick', [],'ZTick', [])
set(gca, 'LooseInset', get(gca, 'TightInset'));
colormap parula
colorbar('off');



% save files.
exportgraphics(fig1,"Solution_for_" + refs + "_Refs.pdf");
exportgraphics(fig2,"Residual_for_" + refs + "_Refs.pdf");