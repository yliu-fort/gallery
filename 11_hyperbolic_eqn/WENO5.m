function [hp, hn] = WENO5(fl, fr, dim)
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
%                               |___________S2__________|
%                               |                       |
%                       |___________S1__________|       |
%                       |                       |       |
%               |___________S0__________|       |       |
%             ..|---o---|---o---|---o---|---o---|---o---|...
%               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
%                                      -|
%                                     i+1/2
%
%
%               |___________S0__________|
%               |                       |
%               |       |___________S1__________|
%               |       |                       |
%               |       |       |___________S2__________|
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
%a=max(abs(dflux(w))); fr=0.5*(flux(w)+a*w); fl=circshift(0.5*(flux(w)-a*w),[0,-1]);
if(nargin<3)
    dim = 1;
end
%% Right Flux
% Choose the positive fluxes, 'v', to compute the left cell boundary flux:
% $u_{i+1/2}^{-}$
frmm = circshift(fr,2,dim);
frm  = circshift(fr,1,dim);
frp  = circshift(fr,-1,dim);
frpp = circshift(fr,-2,dim);

% Polynomials
p0n = (2*frmm - 7*frm + 11*fr)/6;
p1n = ( -frm  + 5*fr  + 2*frp)/6;
p2n = (2*fr   + 5*frp - frpp )/6;

% Smooth Indicators (Beta factors)
B0n = 13/12*(frmm-2*frm+fr  ).^2 + 1/4*(frmm-4*frm+3*fr).^2; 
B1n = 13/12*(frm -2*fr +frp ).^2 + 1/4*(frm-frp).^2;
B2n = 13/12*(fr  -2*frp+frpp).^2 + 1/4*(3*fr-4*frp+frpp).^2;

% Constants
d0n = 1/10; d1n = 6/10; d2n = 3/10; epsilon = 1e-6;

% Alpha weights 
alpha0n = d0n./(epsilon + B0n).^2;
alpha1n = d1n./(epsilon + B1n).^2;
alpha2n = d2n./(epsilon + B2n).^2;
alphasumn = alpha0n + alpha1n + alpha2n;

% ENO stencils weigths
w0n = alpha0n./alphasumn;
w1n = alpha1n./alphasumn;
w2n = alpha2n./alphasumn;

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
hn = w0n.*p0n + w1n.*p1n + w2n.*p2n;

%% Left Flux 
% Choose the negative fluxes, 'u', to compute the left cell boundary flux:
% $u_{i-1/2}^{+}$ 
flmm = circshift(fl,2,dim);
flm  = circshift(fl,1,dim);
flp  = circshift(fl,-1,dim);
flpp = circshift(fl,-2,dim);

% Polynomials
p0p = ( -flmm + 5*flm + 2*fl  )/6;
p1p = ( 2*flm + 5*fl  - flp   )/6;
p2p = (11*fl  - 7*flp + 2*flpp)/6;

% Smooth Indicators (Beta factors)
B0p = 13/12*(flmm-2*flm+fl  ).^2 + 1/4*(flmm-4*flm+3*fl).^2; 
B1p = 13/12*(flm -2*fl +flp ).^2 + 1/4*(flm-flp).^2;
B2p = 13/12*(fl  -2*flp+flpp).^2 + 1/4*(3*fl -4*flp+flpp).^2;

% Constants
d0p = 3/10; d1p = 6/10; d2p = 1/10; epsilon = 1e-6;

% Alpha weights 
alpha0p = d0p./(epsilon + B0p).^2;
alpha1p = d1p./(epsilon + B1p).^2;
alpha2p = d2p./(epsilon + B2p).^2;
alphasump = alpha0p + alpha1p + alpha2p;

% ENO stencils weigths
w0p = alpha0p./alphasump;
w1p = alpha1p./alphasump;
w2p = alpha2p./alphasump;

% Numerical Flux at cell boundary, $u_{i-1/2}^{+}$;
hp = w0p.*p0p + w1p.*p1p + w2p.*p2p;

%hn = circshift(hn,1,dim);
%hp = circshift(hp,1,dim);
% Compute finite volume residual term, df/dx.
%res = (hp-circshift(hp,1,dim)+hn-circshift(hn,1,dim));

% Lax friedrichs flux
%LF = 0.5*(F(hn)+F(hp)-abs(a).*(hp-hn));
end