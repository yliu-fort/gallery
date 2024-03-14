function [hp, hn] = WENO3(v,u, dim)
% *************************************************************************
% Input: u(i) = [u(i-2) u(i-1) u(i) u(i+1) u(i+2)];
% Output: res = df/dx;
%
% Based on:
% C.W. Shu's Lectures notes on: 'ENO and WENO schemes for Hyperbolic
% Conservation Laws' 
%
% coded by Manuel Diaz, 02.10.2012, NTU Taiwan.
% *************************************************************************
%
% Domain cells (I{i}) reference:
%
%                |           |   u(i)    |           |
%                |  u(i-1)   |___________|           |
%                |___________|           |   u(i+1)  |
%                |           |           |___________|
%             ...|-----0-----|-----0-----|-----0-----|...
%                |    i-1    |     i     |    i+1    |
%                |-         +|-         +|-         +|
%              i-3/2       i-1/2       i+1/2       i+3/2
%
% ENO stencils (S{r}) reference:
%
%
%                               |______S1_______|
%                               |               |
%                       |______S0_______|       |
%             ..|---o---|---o---|---o---|---o---|---o---|...
%               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
%                                      -|
%                                     i+1/2
%
%
%                       |______S1_______|
%                       |               |
%                       |       |______S0_______|
%             ..|---o---|---o---|---o---|---o---|---o---|...
%               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
%                               |+
%                             i-1/2
%
% WENO stencil: S{i} = [ I{i-2},...,I{i+2} ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: by using circshift over our domain, we are implicitly creating
% favorable code that includes periodical boundary conditions. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lax-Friedrichs Flux Splitting
%a=max(abs(dflux(w))); v=0.5*(flux(w)+a*w); u=circshift(0.5*(flux(w)-a*w),[0,-1]);
S = 0; % source term
if(nargin<3)
    dim = 1;
end
%% Right Flux
% Choose the positive fluxes, 'v', to compute the left cell boundary flux:
% $u_{i+1/2}^{-}$
vm  = circshift(v,1,dim);
vp  = circshift(v,-1,dim);

% Polynomials
p0n = (-vm + 3*v)/2;
p1n = ( v  + vp )/2;

% Smooth Indicators (Beta factors)
B0n = (vm-v).^2; 
B1n = (v-vp).^2;

% Constants
d0n = 1/3; d1n = 2/3; epsilon = 1E-6;

% Alpha weights 
alpha0n = d0n./(epsilon + B0n).^2;
alpha1n = d1n./(epsilon + B1n).^2;
alphasumn = alpha0n + alpha1n;

% ENO stencils weigths
w0n = alpha0n./alphasumn;
w1n = alpha1n./alphasumn;

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
hn = w0n.*p0n + w1n.*p1n;

%% Left Flux 
% Choose the negative fluxes, 'u', to compute the left cell boundary flux:
% $u_{i-1/2}^{+}$ 
um  = circshift(u,1,dim);
up  = circshift(u,-1,dim);

% Polynomials
p0p = ( um + u )/2;
p1p = (3*u - up)/2;

% Smooth Indicators (Beta factors)
B0p = (um-u).^2; 
B1p = (u-up).^2;

% Constants
d0p = 2/3; d1p = 1/3; epsilon = 1E-6;

% Alpha weights 
alpha0p = d0p./(epsilon + B0p).^2;
alpha1p = d1p./(epsilon + B1p).^2;
alphasump = alpha0p + alpha1p;

% ENO stencils weigths
w0p = alpha0p./alphasump;
w1p = alpha1p./alphasump;

% Numerical Flux at cell boundary, $u_{i-1/2}^{+}$;
hp = w0p.*p0p + w1p.*p1p;
hn = circshift(hn,1,dim);
hp = circshift(hp,-1,dim);

%% Compute finite volume residual term, df/dx.
%res = (hp-circshift(hp,1,dim)+hn-circshift(hn,1,dim)) - S;
%hn = circshift(hn,1,dim);
%hp = circshift(hp,1,dim);
end