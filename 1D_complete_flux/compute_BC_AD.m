function g = compute_BC_AD(x,BCtype)

% impose boundary conditions, given an exact solution
u = exact_solution_AD(x);

if BCtype ==1 %Dirichlet
    g = u;  
end