function [xK,xS,delta_x] = generate_mesh(Nx)

%generate a (uniform) mesh with Nx interior cells
xK = linspace(0,1,Nx+2);
delta_x = xK(2)-xK(1); %since the mesh is uniform, delta_x = xK(j+1)-xK(j)
xS = 0.5*(xK(1:Nx+1) + xK(2:Nx+2)); % interfaces are at midpoints
