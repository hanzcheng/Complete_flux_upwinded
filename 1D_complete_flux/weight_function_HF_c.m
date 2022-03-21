function [w_j,w_j1] = weight_function_HF_c(V,D,delta_x)
%coefficients of c_j and c_{j+1}, homogeneous flux


%Peclet number
Pe = V ./D .* delta_x;


w_j(Pe>0) = -V(Pe>0) ./ (exp(-Pe(Pe>0))-1);
w_j(Pe<=0) = -V(Pe<=0) .* exp(Pe(Pe<=0)) ./ (1-exp(Pe(Pe<=0)));

w_j1(Pe>0) = -V(Pe>0) .* exp(-Pe(Pe>0)) ./ (1-exp(-Pe(Pe>0)));
w_j1(Pe<=0) = -V(Pe<=0) ./ (exp(Pe(Pe<=0))-1);


