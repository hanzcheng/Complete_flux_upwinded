function [c] = solve_drift_diffusion(Nx,D,mu,BC) %homogeneous flux scheme

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


%% compute the velocity 
%at the cell faces
V_face = compute_velocity(mu,xS);

%% coefficients associated with the homogeneous flux, pcwise constant approximation of V
[w_j,w_j1] = weight_function_HF_c(V_face,D,delta_x);
% w_j -> coefficient of c_j
% w_j1 -> coefficient of c_{j+1}


%% Assemble A such that A*c = f_{j+1/2}-f_{j-1/2}
for i=2:Nx+1
    A(i,i-1)= -w_j(i-1);
    A(i,i) =  w_j(i)-w_j1(i-1);
    A(i,i+1)= w_j1(i);
end
%
b(2:Nx+1) = delta_x*compute_source_AD(mu,D,xK(2:Nx+1));

c = A\b;
% end
