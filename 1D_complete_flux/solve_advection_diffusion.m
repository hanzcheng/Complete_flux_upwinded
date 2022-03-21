%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the advection-diffusion equation
% (mu*V*c - Dc')' = s on (0,1), 
% Dirichlet boundary conditions c(0)= c_0, c(1)= c_1
% unknown: c, source term: s
% parameters:
% mu -> mobility coefficient
% V -> velocity
% D -> diffusion

% finite volume method, with unknown c at the nodes at x_j,x_{j+1},
% fluxes are computed at the interface x_{j+1/2}
% o-----|------o
% j   j+1/2   j+1

Nx=20; %Nx = number of cells in the (initial) mesh

% generate a uniform mesh for (0,1) with
[xK,xS,delta_x] = generate_mesh(Nx);
% xK -> nodal points, Nx+2 (Nx interior points + 2 boundary pts) points
% xS -> interface points
% delta_x -> x_{j+1}-x_j

%% parameters
mu = 1;
D = 1e-8;
global epps;
epps = D;
global V_field;
V_field = 3; %3 choices: 1- constant
%           2- linear
%           3- sine function
global testCase;
testCase = 5; %(each test case corresponds to a different exact solution c)
%see exact_solution_AD

%% boundary conditions
BC = [1,1]; %pure Dirichlet BC

%% number of refinements 
nbRefine = 6;
r_param = 2; %refinement parameter, each refinement contains Nx*r_param cells

%initialize variables for measuring error 
err_c = zeros(nbRefine,1);
err_cCF = zeros(nbRefine,1);
err_cCF_up = zeros(nbRefine,1);

%initialize variables for measuring order of convergence, L2 norm
ocl2_c = zeros(nbRefine-1,1);
ocl2_cCF = zeros(nbRefine-1,1);
ocl2_cCF_new = zeros(nbRefine-1,1);


%% solve the drift diffusion equation for nbRefine refinements

for i=1:nbRefine
    %solve system via Homogeneous flux scheme
    c_HF = solve_drift_diffusion(Nx,D,mu,BC);
    %solve system via complete flux, pcwise constant approximation for V
    c_CF = solve_drift_diffusion_CF_std(Nx,D,mu,BC);
    %solve system via complete flux, pcwise linear approximation for V, upwind adjusted
    c_CF_up = solve_drift_diffusion_CF_up(Nx,D,mu,BC);
    
    %% measure errors for convergence test
    [xK] = generate_mesh(Nx);
    c_ex = exact_solution_AD(xK)';
    err_c(i) = norm(c_HF-c_ex)/norm(c_ex);
    err_cCF(i) = norm(c_ex-c_CF)/norm(c_ex);
    err_cCF_up(i) = norm(c_ex - c_CF_up)/norm(c_ex);
    if i>1
        ocl2_c(i) = log(err_c(i-1)/err_c(i))/log(r_param);
        ocl2_cCF(i) = log(err_cCF(i-1)/err_cCF(i))/log(r_param);
        ocl2_cCF_new(i) = log(err_cCF_up(i-1)/err_cCF_up(i))/log(r_param);
    end
    
    %% plot numerical solutions and compare with exact solution
    if mod(i,2)==1
        figure();
        plot(xK,c_ex);
        hold on
        plot(xK,c_HF,'o');
        legend({'exact solution','HF'},'Location','northeast')
        title(['Solution profiles, N=' num2str(Nx)])
        hold off
        
        figure();
        plot(xK,c_ex);
        hold on
        plot(xK,c_CF,'o');
        plot(xK,c_CF_up,'+');
        legend({'exact solution','complete flux','upwinded complete flux'},'Location','northeast')
        title(['Solution profiles, N=' num2str(Nx)])
        hold off
    end
    
    %prepare for next refinement
    Nx = Nx*r_param;
end
