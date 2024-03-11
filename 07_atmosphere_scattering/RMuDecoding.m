function [r, mu] = RMuDecoding(atmosphere, u_r,u_mu)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Rg = atmosphere.Rg;
Rt = atmosphere.Rt;
H = sqrt(Rt^2 - Rg^2);

rho = u_r.*H;
r = sqrt(rho.^2 + Rg.^2);

d_min = Rt - r;
d_max = rho + H;
d = d_min + u_mu .* (d_max - d_min);
mu = (H.^2 - rho.^2 - d.^2)./(2*r.*d);
mu(d == 0) = 1.0;

mu = ClampCosine(mu);
end