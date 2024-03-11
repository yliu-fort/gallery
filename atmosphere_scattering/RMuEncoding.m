function [u_r, u_mu] = RMuEncoding(atmosphere, r,mu)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Rg = atmosphere.Rg;
Rt = atmosphere.Rt;
rho = SafeSqrt(r.^2 - Rg^2);
H = sqrt(Rt^2 - Rg^2);

u_r = rho./H;
d = DistanceToTopAtmosphereBoundary(atmosphere, r, mu);
d_min = Rt - r;
d_max = rho + H;
u_mu = (d - d_min) ./ (d_max - d_min);

end