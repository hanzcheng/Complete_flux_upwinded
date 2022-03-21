function [ue,ue_x,ue_xx] = exact_solution_AD(x_K) 
%exact solution of the advection diffusion equation (and its derivatives)
global testCase;
global epps;
if testCase == 1
ue = ones(size(x_K));
ue_x = zeros(size(x_K));
ue_xx = zeros(size(x_K));
elseif testCase == 2
    ue = x_K;
    ue_x = ones(size(x_K));
    ue_xx = zeros(size(x_K));
elseif testCase == 3   
    ue = sin(pi*x_K);
    ue_x = pi*cos(pi*x_K);
    ue_xx = -pi^2*sin(pi*x_K);
elseif testCase == 4
    a = 1e-1;
    ue = exp( (x_K-1) ./ a );
    ue_x = 1/a * exp( (x_K-1) ./ a);
    ue_xx = 1/a^2 * exp( (x_K-1) ./ a);
elseif testCase == 5
    a = 0.2;
    num = exp((x_K-1)/epps)-exp(-1/epps);
    denum = 1-exp(-1/epps);
    ue = a*sin(pi*x_K) + num ./ denum;
    ue_x = a*pi*cos(pi*x_K) + 1/epps .* exp((x_K-1)/epps) ./ denum;
    ue_xx = -a*pi^2*sin(pi*x_K) + (1/epps^2) .* exp((x_K-1)/epps) ./ denum;
end

    