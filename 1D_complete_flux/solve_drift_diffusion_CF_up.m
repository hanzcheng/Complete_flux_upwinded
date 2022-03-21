function [c] = solve_drift_diffusion_CF_up(Nx,D,mu,BC)
%complete flux scheme, piecewise linear approximation of velocity, upwind adjusted

%% generate mesh
[xK,xS,delta_x] = generate_mesh(Nx);

%% Set boundary conditions
if length(BC)>1
    BCL = BC(1);
    BCR = BC(2);
else
    BCL = BC;
    BCR = BC;
end

%% initialize unknown, matrix system, and source term
c = zeros(Nx+2,1); %Nx finite volumes corresponds to Nx+2 endpoints
A= zeros(Nx+2,Nx+2); %2 boundary points, and Nx interior points
% A is the matrix such that the 1st to Nx+1th entries of
% A*c approximates (Vc'-Dc')
b = zeros(Nx+2,1);

%% impose boundary conditions
%left BC
if BCL==1 %Dirichlet BC
    A(1,1) = 1;
end
b(1) = compute_BC_AD(xK(1),BCL);

%right BC
if BCR==1 %Dirichlet BC
    A(Nx+2,Nx+2)=1;
end
b(Nx+2) = compute_BC_AD(xK(Nx+2),BCR);


%% compute the velocity and its derivative at the cell faces
[V_face,V_dash_face] = compute_velocity(mu,xS);

%% coefficients associated with the homogeneous flux, upwind adjusted V
[w_j,w_j1] = weight_function_HF_up(V_face,V_dash_face,D,delta_x);
% w_j -> coefficient of c_j
% w_j1 -> coefficient of c_{j+1}

%% coefficients associated with the inhomogeneous flux, upwind adjusted V

%compute the Peclet number and the correction term Q
Peclet = V_face/D*delta_x;
Q = V_dash_face/D*delta_x^2/2;
alfa = ones(size(Peclet));
alfa_UB = 0.5*abs(Peclet ./ Q);%upper bound of alpha
alfa = min(alfa_UB,alfa); %factor on [0,1] to ensure that Q does not perform an over-correction
alfa(abs(Peclet)<=10) = 0; %if not in advection-dominated regime, take alpha=0


%compute the modified Peclet number
Peclet_mod = Peclet;
Peclet_mod(Peclet>=0) = Peclet_mod(Peclet>=0)- alfa(Peclet>=0) .* Q(Peclet>=0);
Peclet_mod(Peclet<0) = Peclet_mod(Peclet<0) + alfa(Peclet<0) .* Q(Peclet<0);

Qj = -alfa .* Q;
Qj1 = alfa .* Q;

Qj(Peclet>=0) = 0.25*Qj(Peclet>=0);
Qj1(Peclet>=0) = -0.75*Qj1(Peclet>=0);

Qj(Peclet<0) = -1.25*Qj(Peclet<0);
Qj1(Peclet<0) = -0.25*Qj1(Peclet<0);

%compute coefficients for s_j and s_{j+1}, respectively for the
%inhomogeneous flux
w_i_curr = weight_function_IF_up(-Peclet_mod,Qj);
w_i_next = weight_function_IF_up(Peclet_mod,Qj1);



%% Assemble A such that A*c = f_{j+1/2}-f_{j-1/2}
for i=2:Nx+1
    A(i,i-1)= -w_j(i-1);
    A(i,i) =  w_j(i)-w_j1(i-1);
    A(i,i+1)= w_j1(i);
end
%
b(2:Nx+1) = delta_x*compute_source_AD(mu,D,xK(2:Nx+1));
b_std = b;

for k=3:Nx
    b(k) = b_std(k)+(b_std(k+1)*w_i_next(k)-b_std(k)*w_i_curr(k))-...
        (b_std(k)*w_i_next(k-1)-b_std(k-1)*w_i_curr(k-1));
end
c = A\b;
% end
