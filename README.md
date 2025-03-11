# Numerical solution of optimal control problems using quadratic transport regularization


This repository contains an implementation for solving optimal control problems with optimal transport regularization.
The algorithms are developed in the preprint
[Numerical solution of optimal control problems using quadratic transport regularization](https://arxiv.org/abs/2503.07105).
We minimize $g(u) + W(u)$,
where $g$ is a smooth function and $W$ is the squared quadratic Monge-Kantorovich distance (also known as Wasserstein distance).

The examples 1-6 can be directly executed. The three solution algorithms starting with `solve_` can also be called independently.

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
If you find this helpful, please cite the paper and the source code as
```bibtex
@online{BorchardWachsmuth2025:1,
  author        = {Nicolas Borchard and Gerd Wachsmuth},
  title         = {Numerical solution of optimal control problems using quadratic transport regularization},
  year          = {2025},
  eprint        = {2503.07105},
  eprinttype    = {arXiv}
}
@misc{BorchardWachsmuth2025:2,
  author        = {Nicolas Borchard and Gerd Wachsmuth},
  title         = {Numerical solution of optimal control problems using quadratic transport regularization},
  year          = {2025},
  doi           = {10.5281/zenodo.???}
}
```
