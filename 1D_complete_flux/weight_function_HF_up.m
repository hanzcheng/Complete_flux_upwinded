function [w_j,w_j1] = weight_function_HF_up(V,V_dash,D,delta_x)
%coefficients of c_j and c_{j+1}, upwind adjusted flux'

%Peclet number
Pe = V ./D .* delta_x;
Q = 0.5 * (V_dash ./ D) .* (delta_x .^ 2); 

alfa = ones(size(Pe));
alfa_UB = 2 ./ delta_x .* abs(V ./ V_dash);%upper bound of alpha 
alfa = min(alfa_UB,alfa); %factor on [0,1] to ensure that Q does not perform an over-correction

%upwind adjusted velocity
Vm = V;
Vm (Pe>=0) = V(Pe>=0) - 0.5* alfa(Pe>=0) .* V_dash(Pe>=0) .* delta_x;
Vm (Pe<0) = V(Pe<0) + 0.5* alfa(Pe<0) .* V_dash(Pe<0) .* delta_x;

%upwind adjusted Peclet number
Pe_adj = Pe;
Pe_adj(Pe>=0) = Pe_adj(Pe>=0)- alfa(Pe>=0) .* Q(Pe>=0);
Pe_adj(Pe<0) = Pe_adj(Pe<0)+ alfa(Pe<0) .* Q(Pe<0);

%if Peclet is positive
w_j(Pe>=0) = -Vm(Pe>=0) ./ (exp(-Pe_adj(Pe>=0))-1);
w_j1(Pe>=0) = - exp(-Pe(Pe>=0)) .* w_j(Pe>=0);

%equivalent expression for the weights if Peclet is negative
w_j(Pe<0) = Vm(Pe<0) ./ (exp(Pe_adj(Pe<0))-1);
w_j1(Pe<0) = -w_j(Pe<0);
w_j(Pe<0) = exp(Pe(Pe<0)) .* w_j(Pe<0);