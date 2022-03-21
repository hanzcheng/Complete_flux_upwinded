upwind-adjusted complete flux scheme for an advection-diffusion equation:

Main file: 

solve_advection_diffusion - solves an advection-diffusion problem (with given velocity) via 
1. homogeneous flux scheme
2. complete flux scheme (piecewise constant approximation of velocity)
3. upwind-adjusted complete flux scheme (piecewise linear approximation of velocity)
comparisons are illustrated via plots, and convergence tests


Initialisation: 
1. generate_mesh - generates a uniform mesh with Nx cells
2. compute_velocity - evaluates the velocity and its derivatives at the given point(s) xK
3. exact_solution_AD - manufactures an exact solution and its derivatives (until 2nd order)
	3a. compute_BC_AD - impose the boundary conditions according to the exact solution
	3b. compute_source_AD - piecewisec constant approximation of the source term, based on the exact solution

Solving: 

1. solve_drift_diffusion - solves the advection-diffusion equation via homogeneous flux (Scharfetter Gummel scheme)
2. solve_drift_diffusion_CF_std - solves the advection-diffusion equation via complete flux scheme (piecewise constant approximation of velocity)
3. solve_drift_diffusion_CF_up - solves the advection-diffusion equation via upwind-adjusted complete flux scheme (piecewise linear approximation of velocity)


Matrix assembly:

1. weight_function_HF_c - coefficients associated to the homogeneous flux (piecewise constant approximation of velocity)
2. weight_function_HF_up - coefficients associated to the homogeneous flux (piecewise linear approximation of velocity)
3. weight_function_IF_c - coefficients associated to the inhomogeneous flux (piecewise constant approximation of velocity)
4. weight_function_IF_up - coefficients associated to the inhomogeneous flux (piecewise linear approximation of velocity)