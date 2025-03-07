%plots two different solutions for non-convex case.

% yd_fun = @(x)(cos(pi/2*x(1,:)).*cos(pi/2*x(2,:)).*exp(x(2,:)) / 2);
% yd2_fun = @(x)(0.6*cos(pi/2*x(1,:)+pi).*sin(pi*x(2,:)) / 2);

clear;


%% prepare situation.
%non-convex problem.
problem.convex = false;
% The desired states
yd_fun = @(x)(cos(pi/2*x(1,:)).*cos(pi/2*x(2,:)).*exp(x(2,:)) / 2);
yd2_fun = @(x)(0.5*cos(pi/2*x(1,:)+pi).*sin(pi*x(2,:)) / 2);
% Refine the mesh
refs=4;

vertices_D = ...
  [-1,-1;
    1,-1;
    1, 1;
   -1, 1];
problem.D.vertices = vertices_D;
vertices_Omega = vertices_D;
problem.F.vertices = vertices_D;
% Weight parameters
problem.alpha = 1e-3;

problem = setup(refs,problem);

ydi = problem.ydi;
yd2i = problem.yd2i;
Mi = problem.Mi;
Ki = problem.Ki;
I4 = problem.I4;



%% give start values for different solutions.
norm1 = .5*((yd2i - ydi)'*(Mi*(yd2i - ydi)));
norm2 = .5*((ydi - yd2i)'*(Mi*(ydi - yd2i)));

p = Ki\(Mi*( norm2 * ( yd2i - ydi) ));

w1 = problem.A + I4*p/problem.alpha;
w2 = problem.A - I4*p/problem.alpha;



%% calculate solutions.
% algo = 1; %Newton
% algo = 2; %Fixed-point
algo = 3; %Fixed-point fixed-step

if algo == 1
		[solution1] = solve_newton(problem,w1);
		[solution2] = solve_newton(problem,w2);
elseif algo == 2
		[solution1] = solve_fixed_point(problem,w1);
		[solution2] = solve_fixed_point(problem,w2);
else
		para.ls_alpha = 3e-2;
		% para.ls_alpha = 7e-2; %no convergence
		para.abstol = 1e1;
		para.maxiter = 30;

		[solution1] = solve_fixed_point_fixed_step(problem,w1,para);
		[solution2] = solve_fixed_point_fixed_step(problem,w2,para);

		para.abstol = 1e-12;

		[solution1] = solve_newton(problem,solution1.w,para);
		[solution2] = solve_newton(problem,solution2.w,para);
end
u1 = solution1.u;
u2 = solution2.u;


%% show and save.
g = problem.g;

fig1 = figure(1); clf
pdeplot(g.p, g.e, g.t, 'xydata', u1, 'zdata', u1, 'mesh', 'on');
% title('control $\bar u_1$', 'interpreter', 'latex','FontSize',32)
set(gca, 'LooseInset', get(gca, 'TightInset'));
xlabel('$x_1$','interpreter', 'latex','FontSize',80,'Position', [mean(xlim) min(ylim) min(u1)]);
ylabel('$x_2$','interpreter', 'latex','FontSize',80,'Position', [min(xlim) mean(ylim) min(u1)]);
set(gca,'XTick',[], 'YTick', [],'ZTick', [])
colormap(parula);
colorbar('off');

fig2 = figure(2); clf
pdeplot(g.p, g.e, g.t, 'xydata', u2, 'zdata', u2, 'mesh', 'on');
% title('control $\bar u_2$', 'interpreter', 'latex','FontSize',32)
set(gca, 'LooseInset', get(gca, 'TightInset'));
xlabel('$x_1$','interpreter', 'latex','FontSize',80,'Position', [mean(xlim) min(ylim) min(u2)]);
ylabel('$x_2$','interpreter', 'latex','FontSize',80,'Position', [min(xlim) mean(ylim) min(u2)]);
set(gca,'XTick',[], 'YTick', [],'ZTick', [])
colormap(parula);
colorbar('off');


Y = problem.I2*(solution1.y - problem.ydi);
fig3 = figure(3); clf
pdeplot(g.p, g.e, g.t, 'xydata', Y, 'zdata', Y, 'mesh', 'on');
% title('residual $y_1 - y_{d,i}$','interpreter','latex','FontSize',32)
xlabel('$x_1$','interpreter', 'latex','FontSize',80,'Position', [mean(xlim) min(ylim) min(Y)]);
ylabel('$x_2$','interpreter', 'latex','FontSize',80,'Position', [min(xlim) mean(ylim) min(Y)]);
set(gca,'XTick',[], 'YTick', [],'ZTick', [])
colormap(parula);
colorbar('off');

Y = problem.I2*(solution2.y - problem.yd2i);
fig4 = figure(4); clf
pdeplot(g.p, g.e, g.t, 'xydata', Y, 'zdata', Y, 'mesh', 'on');
% title('residual $y_2 - y_{d,2,i}$','interpreter','latex','FontSize',32)
xlabel('$x_1$','interpreter', 'latex','FontSize',80,'Position', [mean(xlim) min(ylim) min(Y)]);
ylabel('$x_2$','interpreter', 'latex','FontSize',80,'Position', [min(xlim) mean(ylim) min(Y)]);
set(gca,'XTick',[], 'YTick', [],'ZTick', [])
colormap(parula);
colorbar('off');



exportgraphics(fig1,'Nonconvex_u1.pdf');
exportgraphics(fig2,'Nonconvex_u2.pdf');
exportgraphics(fig3,'Nonconvex_res1.pdf');
exportgraphics(fig4,'Nonconvex_res2.pdf');
