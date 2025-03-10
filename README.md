# Numerical solution of optimal control problems using quadratic transport regularization


This repository contains an implementation of an example for optimal control problems using optimal transport regularization made in the preprint
[Numerical solution of optimal control problems using quadratic transport regularization](arxiv todo).
There $g(u)+W(u)$ is minimized,
where $g$ is a smooth function and $W$ is the squared quadratic Monge-Kantorovich distance (also known as Wasserstein distance).

The examples 1-6 and additionally the three algorithms starting with "solve" can be executed.
Here, we shortly write down what they do.

`solve_fixed_point.m`: Armijo line search algorithm using $r(w)$ as search direction. The equation $r(w) = 0$ is the optimality condition of the problem. It can also be seen as some gradient method. This algorithm is rather slow.

`solve_fixed_point_fixed_step.m`: $r(w)$ is used as search direction together with constant step sizes.

`solve_newton.m`: semismooth Newton method is applied to $r(w) = 0$. This algorithm is quite fast.

`ex_1_optimal_solution.m`: calculates the optimal solution for some instance of the problem class.

`ex_2_transport.m`: shows the transport map related to the optimal solution.

`ex_3_speed_test.m`: compares the speed of the fixed point and Newton method for different levels of refinements of a coarse mesh.

`ex_4_different_refinements.m`: plots a diagram on how the infinity norm of $r(w)$ decreases by every iteration vor different levels of refinement. 

`ex_5_subset_control.m`: the control $u$ only lives on a subset of the mesh.

`ex_6_nonconvex.m`: shows two different local solutions for a nonconvex smooth function $g$.

## Citation
If you find this helpful, please cite our manuscript
```bibtex
@online{BorchardWachsmuth2025:1,
  author        = {},
  title         = {},
  year          = {},
  eprint        = {},
  eprinttype    = {arXiv}
}
```
