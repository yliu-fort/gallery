clear;clc;close all;

% nd reconstruction
A = ndgrid(rand(16,1),rand(64,1),rand(16,1),rand(8,1));
sA = size(A);
A0 = A(:);
A0 = reshape(A0, sA);
assert(all(A(:) == A0(:)));

clear
atmosphere.Rg = 1.0;
atmosphere.Rt = 1.011;

% Debug r, mu encoding/decoding
u_rDist = linspace(0,1,64);
u_muDist = linspace(0,1,256);
[u_r,u_mu] = ndgrid(u_rDist, u_muDist);
[r,mu] = RMuDecoding(atmosphere, u_r, u_mu);
[u_r1,u_mu1] = RMuEncoding(atmosphere, r, mu);

fprintf("r, mu reconstruction error = \n %e, %e.\n", ...
    norm(u_r(:) - u_r1(:)), ...
    norm(u_mu(:) - u_mu1(:)));



u_rDist = linspace(0,1,16);
u_muDist = linspace(0,1,64);
u_musDist = linspace(0,1,16);
u_nuDist = linspace(0,1,4);
[u_r,u_mu,u_mu_s,u_nu] = ndgrid(u_rDist,u_muDist,u_musDist,u_nuDist);

clear
atmosphere.Rg = 1.0;
atmosphere.Rt = 1.011;
sunDir = [cosd(90),0,sind(90)];
viewPos = [0,0,1.0001];
theta = linspace(0, 2*pi, 1);
phi = linspace(pi/2, -pi/2*0.9999, 512);
[T, P]=meshgrid(theta, phi);
I = [];

viewDir = [];
viewDir(:,1) = cos(T(:)).*cos(P(:));
viewDir(:,2) = sin(T(:)).*cos(P(:));
viewDir(:,3) = sin(P(:));

% Process indata
r = norm(viewPos) + 0.*T(:);
mu = dot(viewPos+ 0*T(:),viewDir,2)./vecnorm(viewPos,2,2);
mu_s = dot(viewPos+ 0*T(:), sunDir + 0*T(:),2)./vecnorm(viewPos,2,2);
nu = dot(viewDir, sunDir + 0*T(:),2);


[u_r, u_mu, u_mu_s, u_nu] = RMuMusNuEncoding(atmosphere, r ...
    , mu, mu_s, nu);
[r1, mu1, mu_s1, nu1, ~] = RMuMusNuDecoding(atmosphere, u_r ...
    , u_mu, u_mu_s, u_nu);


fprintf("\nr, mu, mu_s, nu \n reconstruction error = \n %e, \n %e, \n%e, \n%e. \n", ...
    norm(r(:) - r1(:)), ...
    norm(mu(:) - mu1(:)),...
    norm(mu_s(:) - mu_s1(:)), ...
    norm(nu(:) - nu1(:)));