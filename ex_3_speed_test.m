clear;

% Fix random seed
rng(1);

minnrofref = 4;
maxnrofref = 8;

%define basic setting of problems.
basic_problem.D.vertices = ...
            [-1,-1;
              1,-1;
              1, 1;
             -1, 1];
basic_problem.vertices_Omega = basic_problem.D.vertices;
basic_problem.F.vertices = basic_problem.D.vertices;

for i = minnrofref:maxnrofref
	%make problem for i refinements.
	problem_i = setup(i,basic_problem);

	%calculate start values.
	w = problem_i.A + randn(size(problem_i.a,1),1)*1e-7;

	%calculate solutions for different algorithms.
	s1 = solve_fixed_point(problem_i,w);
	s2 = solve_newton(problem_i,w);

	index = 1+i-minnrofref;
	sol_fixed(index)= s1;
	sol_newton(index)= s2;
end

for refs = minnrofref:maxnrofref
	index = 1+refs-minnrofref;

	s1 = sol_fixed(index);
	s2 = sol_newton(index);

	fprintf('%4d & %10d & %10d & %10.2f & %6.3e & %10d & %10.2f & %6.3e \\\\\n', ...
		refs, length(s1.w), ...
		s1.iter, s1.time, s1.time / s1.iter, ...
		s2.iter, s2.time, s2.time / s2.iter);
end

