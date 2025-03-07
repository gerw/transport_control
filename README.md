This program code is an implementation of an example for optimal control problems using optimal transport regularization made in the following paper:
[hier einfügen]
There g(u)+W(u) is minimized,
where g is a smooth function and W is the Monge-Kantorovich distance (also known as Wasserstein distance).

The examples 1-6 and additionally the three algorithms starting with "solve" can be executed.
Here, we shortly write down what they do.

solve_fixed_point: Armijo line search algorithm using r(w) as search direction. The equation r(w) = 0 is the optimality condition of the problem. It can also be seen as some gradient method. This algorithm is rather slow.

solve_fixed_point_fixed_step: r(w) is used as search direction together with constant stepsizes.

solve_newton: semismooth Newton method is applied to r(w). This algorithm is quite fast.

ex_1_optimal_solution: calculates the optimal solution for some instance of the problem class.

ex_2_transport: shows the transport map related to the optimal solution.

ex_3_speed_test: compares the speed of the fixed point and Newton method for different levels of refinements of a coarse mesh.

ex_4_different_refinements: plots a diagram on how the infinity norm of r(w) decreases by every iteration vor different levels of refinement. 

ex_5_subset_control: the control u only lives on a subset of the mesh.

ex_6_nonconvex: shows two different local solutions for a nonconvex smooth function g.
