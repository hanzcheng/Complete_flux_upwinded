function s = compute_source_AD(mu,D,x_K)

%source term of the advection-diffusion equation, given exact solution ue
[ue,ue_x,ue_xx] = exact_solution_AD(x_K);
[V,V_dash] = compute_velocity(mu,x_K);

s =ue .* V_dash  + V .* ue_x - D .* ue_xx;

end