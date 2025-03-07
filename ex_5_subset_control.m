clear;


%define problem structure.
refs = 5;

problem.D.vertices = ...
            [-1,-1;
              1,-1;
              1, 1;
             -1, 1];
problem.vertices_Omega = problem.D.vertices;
problem.F.vertices = ...
  [ 0, 0;
    1, 0;
    1, 1;
    0, 1];

problem.alpha = 1e-3;

problem  = setup(refs,problem);


%solve problem.
solution = solve_newton(problem);


%% plot
u = solution.u;
g = problem.g;
fig = figure(1); clf
pdeplot(g.p, g.e, g.t, 'xydata', u, 'zdata', u, 'mesh', 'on');
xlabel('$x_1$','interpreter', 'latex','FontSize',48,'Position', [mean(xlim) min(ylim) min(u)]);
ylabel('$x_2$','interpreter', 'latex','FontSize',48,'Position', [min(xlim) mean(ylim) min(u)]);
set(gca,'XTick',[], 'YTick', [],'ZTick', [])
colorbar('off')
colormap(parula);
set(gca, 'LooseInset', get(gca, 'TightInset'));
% title('control $\bar u$','interpreter','latex','FontSize',28)

%save file.
exportgraphics(fig,"Controlled_only_on_subset.pdf");

