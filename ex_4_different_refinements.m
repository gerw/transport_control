%%helper script for plotting the value of the infinity norm of r,
%%the infinity norm of the gradient of J as well as J over time for
%%different numbers of refinements of the mesh.
%%Creates files with that data.

clear;

%maximum number of refinements.
minnrofref = 4;
maxnrofref = 9;
solutions = [];
for i = minnrofref:maxnrofref
	problem_i = setup(i);
	%As the discretization changes the dimension of the problem changes.
	%Thus, initial values are different by discretization number.
	w = problem_i.A + randn(size(problem_i.a,1),1)*1e-7;
	solution = solve_newton(problem_i,w);
	solutions = [solutions, solution];
end

%% plot graphs now.
fig = figure(3); clf
hold on;
for i = minnrofref:maxnrofref
	ind = i-minnrofref+1;
	sol_i = solutions(ind);
	% plot(1:sol_i.iter,sol_i.list_dJ);
	plot(1:sol_i.iter,sol_i.list_dJ, 'DisplayName',sprintf('%d',i),'LineWidth',3);
	% plot(1:sol_i.iter,sol_i.list_dJ,'Color',cmap(i,:),'DisplayName');
	text(sol_i.iter,sol_i.list_dJ(end),num2str(i),FontSize=44);
end
xlabel('Iterations $k$', 'interpreter', 'latex','Fontsize',44);
ylabel('$|\nabla J(w_k)|_\infty$', 'interpreter', 'latex','Fontsize',44);
set(gca, 'YScale', 'log',FontSize=44);
% title('$|\nabla J|_\infty$ by iterations and discretization (Newton)', 'interpreter', 'latex');
set(gca, 'LooseInset', get(gca, 'TightInset'));
hold off;
legend show;

%save file.
exportgraphics(fig,"dJinf_" + maxnrofref + "Ref.pdf");


%next graph
fig = figure(4); clf
hold on;
for i = minnrofref:maxnrofref
	ind = i-minnrofref+1;
	sol_i = solutions(ind);
	plot(1:sol_i.iter,sol_i.list_J, 'DisplayName',sprintf('%d',i),'LineWidth',3);
	text(sol_i.iter,sol_i.list_J(end),num2str(i),FontSize=32);
end
xlabel('Iterations $k$', 'interpreter', 'latex','Fontsize',44);
ylabel('$J(w_k)$', 'interpreter', 'latex','Fontsize',44);
set(gca, 'YScale', 'log',FontSize=44);
% title('$J(w_k)$ by iterations and discretization (Newton)', 'interpreter', 'latex');
set(gca, 'LooseInset', get(gca, 'TightInset'));
hold off;
legend show;

%save file.
exportgraphics(fig,"J_" + maxnrofref + "Ref.pdf");


%next graph
fig = figure(5); clf
hold on;
for i = minnrofref:maxnrofref
	ind = i-minnrofref+1;
	sol_i = solutions(ind);
	plot(1:sol_i.iter,sol_i.list_r,'DisplayName',sprintf('%d',i),'LineWidth',3);
	if i == 7
		text(sol_i.iter,sol_i.list_r(end)+2e-7,num2str(i),FontSize=32);
	elseif i == 8
		text(sol_i.iter,sol_i.list_r(end)-3e-8,num2str(i),FontSize=32);
	else
		text(sol_i.iter,sol_i.list_r(end),num2str(i),FontSize=32);
	end
end
xlabel('Iterations $k$', 'interpreter', 'latex');
ylabel('$|r(w_k)|_\infty$', 'interpreter', 'latex');
set(gca, 'YScale', 'log',FontSize=44);
% title('$|r(w_k)|_\infty$ by iterations and discretization (Newton)', 'interpreter', 'latex');
set(gca, 'LooseInset', get(gca, 'TightInset'));
hold off;
legend show;

%save file.
exportgraphics(fig,"rinf_" + maxnrofref + "Ref.pdf");
