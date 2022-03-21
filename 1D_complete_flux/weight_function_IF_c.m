function w_src = weight_function_IF_c(P_j)
%coefficients needed for inhomogeneous flux 

tol=1e-6;

Pjn = find(P_j<0);
num_src(Pjn) = exp(0.5*P_j(Pjn))-0.5*(P_j(Pjn)) -1;
denum_src(Pjn) = (P_j(Pjn)) .* (exp(P_j(Pjn))-1); 

Pjp = find(P_j>0);
num_src(Pjp) = exp(-P_j(Pjp)*0.5)- exp(-P_j(Pjp)) .* (0.5*(P_j(Pjp))+1);
denum_src(Pjp) = (P_j(Pjp)) .* (1-exp(-P_j(Pjp))); 

w_src = num_src ./ denum_src;

% if Peclet is close to 0, use the limit as P->0 for W
w_src(abs(P_j)<tol) = 0.125;
