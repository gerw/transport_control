function [] = do_plots(problem,solution)
% Plots the solution

g = problem.g;
u = solution.u;
ydi = problem.ydi;
Ki = problem.Ki;
Mi = problem.Mi;


if problem.n > 1e5
	plotmesh = 'off';
else
	plotmesh = 'on';
end


if problem.plot_u
	%u is the height of the Diracs on a_i. We plot it by filling all other
	%values affinely on triangles. We also lumped it.
	figure(1); clf
	repu = u ./ sum(problem.M,2);
	% The next line sets the control to zero on the boundary.
	% repu = problem.I2*problem.I2'*repu;
	pdeplot(g.p, g.e, g.t, 'xydata', repu, 'zdata', repu, 'mesh', plotmesh);
	% pdeplot(g.p, g.e, g.t, 'xydata', repu, 'zdata', repu);
	% pdeplot(g.p, g.e, g.t, 'xydata', log(u), 'zdata', log(u), 'mesh', plotmesh);
	title('control $\bar u$', 'interpreter', 'latex')
	xlabel('$x_1$','interpreter', 'latex','FontSize',32,'Position', [mean(xlim) min(ylim) min(u)]);
	ylabel('$x_2$','interpreter', 'latex','FontSize',32,'Position', [min(xlim) mean(ylim) min(u)]);
	set(gca,'XTick',[], 'YTick', []);
	set(gca, 'Fontsize',28);
	colormap parula;
end



y = Ki\( problem.I2'*u );
Y = problem.I2*(ydi-y);

if problem.plot_residuum
	figure(2); clf
	pdeplot(g.p, g.e, g.t, 'xydata', Y, 'zdata', Y, 'mesh', plotmesh);
	title('residual')
	xlabel('$x_1$','interpreter', 'latex','FontSize',32,'Position', [mean(xlim) min(ylim) min(Y)]);
	ylabel('$x_2$','interpreter', 'latex','FontSize',32,'Position', [min(xlim) mean(ylim) min(Y)]);
	set(gca,'XTick',[], 'YTick', []);
	set(gca, 'Fontsize',28);
	colormap parula;
end

if problem.plot_p
	figure(3); clf
	if true
		% Plot adjoint
		P = problem.I2*(Ki\(Mi*(y - ydi)));


		figure(3); clf
		pdeplot(g.p, g.e, g.t, 'xydata', P, 'zdata', P, 'mesh', plotmesh);
		title('adjoint $\bar p$', 'interpreter', 'latex')
		xlabel('$x_1$','interpreter', 'latex','FontSize',32,'Position', [mean(xlim) min(ylim) min(P)]);
		ylabel('$x_2$','interpreter', 'latex','FontSize',32,'Position', [min(xlim) mean(ylim) min(P)]);
		set(gca,'XTick',[], 'YTick', []);
		set(gca, 'Fontsize',28);
		colormap parula;
	else
		% Plot modified adjoint
		p = problem.I2*(Ki\(Mi*(y - ydi)));



		P = p/problem.alpha + sum(g.p.^2)'/2;
		pdeplot(g.p, g.e, g.t, 'xydata', P, 'zdata', P, 'mesh', plotmesh);
		title('scaled adjoint plus norm squared, i.e. $A+p/\alpha = w$', 'interpreter', 'latex')
		xlabel('$x_1$','interpreter', 'latex','FontSize',32,'Position', [mean(xlim) min(ylim) min(P)]);
		ylabel('$x_2$','interpreter', 'latex','FontSize',32,'Position', [min(xlim) mean(ylim) min(P)]);
		set(gca,'XTick',[], 'YTick', []);
		set(gca, 'Fontsize',28);
		colormap parula;
	end
end

%Additionally in the non-convex case
if problem.plot_y
	Y = problem.I2*( problem.Ki\( problem.I2' * u ) );
	figure(4); clf
	pdeplot(g.p, g.e, g.t, 'xydata', Y, 'zdata', Y, 'mesh', plotmesh);
	title('$y$','interpreter', 'latex')
	xlabel('$x_1$','interpreter', 'latex','FontSize',32,'Position', [mean(xlim) min(ylim) min(Y)]);
	ylabel('$x_2$','interpreter', 'latex','FontSize',32,'Position', [min(xlim) mean(ylim) min(Y)]);
	set(gca,'XTick',[], 'YTick', []);
	set(gca, 'Fontsize',28);
	colormap parula;
end

if problem.plot_yd
	Ydi = problem.I2*problem.ydi;
	figure(5); clf
	pdeplot(g.p, g.e, g.t, 'xydata', Ydi, 'zdata', Ydi, 'mesh', plotmesh);
	title('$y_{d,i}$','interpreter', 'latex')
	xlabel('$x_1$','interpreter', 'latex','FontSize',32,'Position', [mean(xlim) min(ylim) min(Ydi)]);
	ylabel('$x_2$','interpreter', 'latex','FontSize',32,'Position', [min(xlim) mean(ylim) min(Ydi)]);
	set(gca,'XTick',[], 'YTick', []);
	set(gca, 'Fontsize',28);
	colormap parula;
end

if problem.plot_yd2
	Yd2i = problem.I2*problem.yd2i;
	figure(6); clf
	pdeplot(g.p, g.e, g.t, 'xydata', Yd2i, 'zdata', Yd2i, 'mesh', plotmesh);
	title('$y_{d,2,i}$','interpreter', 'latex')
	xlabel('$x_1$','interpreter', 'latex','FontSize',32,'Position', [mean(xlim) min(ylim) min(Yd2i)]);
	ylabel('$x_2$','interpreter', 'latex','FontSize',32,'Position', [min(xlim) mean(ylim) min(Yd2i)]);
	set(gca,'XTick',[], 'YTick', []);
	set(gca, 'Fontsize',28);
	colormap parula;
end


drawnow;
end
