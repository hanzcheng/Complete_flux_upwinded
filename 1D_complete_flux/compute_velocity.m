function [V,V_dash] = compute_velocity(mu,xK)
% velocity = mu*V
% compute the velocity and its derivative at the point(s) xK

global V_field

if V_field == 1 %constant velocity field 
    V = ones(size(xK));
    V_dash = zeros(size(xK));
elseif V_field == 2 %linear velocity field (V<0)
    
    V = (xK-1.5);
    V_dash = ones(size(xK));
    
elseif V_field == 3 %sine function (V>0)
    V = (1+0.95*sin(pi*xK));
    V_dash = 0.95*pi*cos(pi*xK);
end

V = mu*V;
V_dash = mu*V_dash;