clear; clc; close all;

%% Inviscid Burgers
gamma = 1.4;
t0 = 0;
tf = 0.2;
idim = 100;
ti = 0.002;
flux = zeros(idim+2,3);

U = zeros(idim+2,3);

dx = 1.0/idim;

U = test1(U,idim+2, gamma);

%% Main program
t = t0;
nti = 1;
tif = ti*nti;

sol.rho = [];
sol.u = [];
sol.p = [];
sol.e = [];
sol.m = [];

for i = 1:100
    U = BoundaryCorrection(U);
    dt = cflcon(0.4, U(:,2)./U(:,1), getSoundSpeed(U, gamma), dx, t, tif);
    t = t + dt;
    %flux_op = @fluxes;
    %flux_op = @KT_fluxes;
    flux_op = @WENO_fluxes;
    %flux_op = @WENO_fluxes_fv;
    %flux_op = @MUSCL_fluxes;
    %flux_op = @Rusanov_fluxes;
    %flux_op = @AUSMDV_fluxes;
    %U = update(U, flux, dt, dx);
    %U =  U - dt/dx*flux;


    Uo = U;

    flux = flux_op(flux,U,dx,dt,gamma,3);
    U =  Uo - dt/dx*flux;
    U = BoundaryCorrection(U);

    flux = flux_op(flux,U,dx,dt,gamma,3);
    U = 0.75*Uo+0.25*(Uo-dt/dx*flux);
    U = BoundaryCorrection(U);

    flux = flux_op(flux,U,dx,dt,gamma,3);
    U = (Uo+2*(U-dt/dx*flux))/3;
    U = BoundaryCorrection(U);

    
    %if(t == tif)
        sol.rho = [sol.rho U(:,1)];
        sol.u = [sol.u U(:,2)./U(:,1)];
        sol.m = [sol.m U(:,2)./U(:,1)./getSoundSpeed(U, gamma)];
        sol.p = [sol.p getPressure(U,gamma)];
        sol.e = [sol.e getInternalEnergy(U)];
        nti = nti+1;
        tif = ti*nti;
        
        disp(t)
    %end

end

%% Plot
figure(1)
field = 'rho';
surf(sol.(field).','EdgeColor','none');
view(2), axis tight;
colormap gray;

%% Plot 2
figure(2)
subplot(2,2,1)
plot(sol.('rho')(:,end),'.-'),hold on,title('Rho')
subplot(2,2,2)
plot(sol.('m')(:,end),'.-'),hold on,title('mach')
subplot(2,2,3)
plot(sol.('p')(:,end),'.-'),hold on,title('p')
subplot(2,2,4)
plot(sol.('e')(:,end),'.-'),hold on,title('e')

%% MUSCL_Hancock
function phi = mc(a,b)
d = (a > 0) & (b > 0) -(a < 0) & (b < 0) ;
phi = d.*min(min(2*abs(a),2*abs(b)),abs(a+b)/2);
end

function flux = WENO_fluxes(flux,U,dx,dt,gamma,nvars)
% CFL < 0.5 to eliminate strips
% Data reconstruction X
W = ToPrimitiveVars(U,gamma);
a = max(abs(W(:,2))+abs(getSoundSpeed(U, gamma)));

%sflux=FV_primWise_WENO5LF1d(a,U.',2);
%sflux = sflux.';
%a = max(max(max(abs(W(:,:,2))+abs(getSoundSpeed(U, gamma)))));

[HWL,HWR] = WENO5(circshift(W,-1,1),W, 1);
HUL = ToConservedVars(HWL,gamma);
HUR = ToConservedVars(HWR,gamma);
F = 0.5*(compute_flux(HWR,gamma)+compute_flux(HWL,gamma) - abs(a).*(HUL-HUR));
flux = F - circshift(F,1,1);

end

function flux = KT_fluxes(flux,U,dx,dt,gamma,nvars)
% CFL < 0.5 to eliminate strips
% Data reconstruction X
W = ToPrimitiveVars(U,gamma);
Wm1 = circshift(W, 1,1); %i-1
Wp1 = circshift(W,-1,1); %i+1

dW = mc(Wp1 - W, W - Wm1);

WL = W + 0.5*dW;
WL = circshift(WL,1,1);
WR = W - 0.5*dW;

UL = ToConservedVars(WL,gamma);
UR = ToConservedVars(WR,gamma);

AL = getSoundSpeed(UL, gamma);
AR = getSoundSpeed(UR, gamma);
a = max(abs(WL(:,2))+AL,abs(WR(:,2))+AR);
flux = 0.5*(compute_flux(WR,gamma)+compute_flux(WL,gamma) - a.*(UR-UL));
flux = flux(2:end,:);
flux = diff(flux);
flux = [zeros(1,3);flux;zeros(1,3)];

end

function flux = MUSCL_fluxes(flux,U,dx,dt,gamma,nvars)
% Data reconstruction
%dU = diff(U);
%dl = [zeros(1,nvars);dU(1:end-1,:);zeros(1,nvars)];
%dr = [zeros(1,nvars);dU(2:end  ,:);zeros(1,nvars)];

% Minmod limiter
%beta = 1;
%di = (dr > 0).*max(0, max(min(beta*dl,dr), min(dl,beta*dr))) ...
%   + (dr < 0).*min(0, min(max(beta*dl,dr), max(dl,beta*dr)));

W = ToPrimitiveVars(U,gamma);
Wm1 = circshift(W, 1,1); %i-1
Wp1 = circshift(W,-1,1); %i+1


WL = W + 0.5*(W - Wm1).*max(0,min(1,(Wp1-W)./(W - Wm1)));
WL = circshift(WL,1,1);
WR = W + 0.5*(W - Wp1).*max(0,min(1,(Wm1-W)./(W - Wp1)));

% Evolution
%UL = U - 0.5.*di;
%UR = U + 0.5.*di;

% fluxL = compute_flux(ToPrimitiveVars(UL, gamma),gamma);
% fluxR = compute_flux(ToPrimitiveVars(UR, gamma),gamma);
% 
% WL = ToPrimitiveVars(UL + 0.5*dt/dx.*(fluxL - fluxR),gamma);
% WR = ToPrimitiveVars(UR + 0.5*dt/dx.*(fluxL - fluxR),gamma);

fluxL = compute_flux(WL,gamma);
fluxR = compute_flux(WR,gamma);

WL = WL - 0.5*dt/dx.*(fluxR - fluxL);
WR = WR - 0.5*dt/dx.*(fluxR - fluxL);

% Solve Riemann problem
% Call riemann solver for inviscid burger

riemann_solver_t = @(rhol,ul,pl,rhor,ur,pr)(riemann_solver(rhol,ul,pl,rhor,ur,pr,gamma));
[~, ~, sampler] = arrayfun(riemann_solver_t,...
    WR(1:end-1,1),WR(1:end-1,2),WR(1:end-1,3),...
    WL(2:end,1),WL(2:end,2),WL(2:end,3),...
    'UniformOutput',false);

% Compute Godunov intercell flux
for i = 1:numel(sampler)
    [rs,us,ps] = feval(sampler{i},0,dt);
    flux(i,:) = compute_flux([rs,us,ps],gamma);
end

%{
% Roe Solver
riemann_solver_t = @(rhol,ul,pl,rhor,ur,pr)(roe_solver(rhol,ul,pl,rhor,ur,pr,gamma));
[flux(2:end,1),flux(2:end,2),flux(2:end,3)] = arrayfun(riemann_solver_t,...
    WR(1:end-1,1),WR(1:end-1,2),WR(1:end-1,3),...
    WL(2:end,1),WL(2:end,2),WL(2:end,3));
%}
flux = flux(2:end,:);
flux = diff(flux);
flux = [zeros(1,3);flux;zeros(1,3)];
mustBeReal(flux);

end

%% Helper functions
function umag = getVelMagSqr(U)
    umag = sum((U(:,2:(end-1))./U(:,1)).^2,2);
end

function e = getInternalEnergy(U)
    e = (U(:,end) - 0.5*U(:,1).*getVelMagSqr(U))./U(:,1);
end

function p = getPressure(U, gamma)
    p = (gamma-1)*(U(:,end) - 0.5*U(:,1).*getVelMagSqr(U));
end

function a = getSoundSpeed(U, gamma)
    %a = sqrt(gamma*(gamma-1)*getInternalEnergy(U));
    a = sqrt(gamma*getPressure(U,gamma)./U(:,1));
end

function W = ToPrimitiveVars(U, gamma)
W(:,1) = U(:,1);
W(:,2) = U(:,2)./U(:,1);
W(:,3) = (U(:,3) - 0.5.*U(:,2).^2./U(:,1))*(gamma-1);
end

function U = ToConservedVars(W, gamma)
U(:,1) = W(:,1);
U(:,2) = W(:,1).*W(:,2);
U(:,3) = 0.5*W(:,1).*W(:,2).^2 + W(:,3)/(gamma-1);
end

function U = test1(U,idim, gamma)
rhol = 1.0;ul = 0.75;pl = 1.0;
rhor = 0.125;ur =  0.0;pr = 0.1;

xpos = -.5;
for i = 1: idim
    xpos = xpos + 1.0/idim;
    
    if(xpos < 0)
        U(i,1) = rhol;
        U(i,2) = rhol*ul;
        U(i,3) = rhol*(0.5*ul^2+pl/(gamma-1)/rhol);
    else
        U(i,1) = rhor;
        U(i,2) = rhor*ur;
        U(i,3) = rhor*(0.5*ur^2+pr/(gamma-1)/rhor);
    end
end

end

% two shock
function U = test2(U,idim, gamma)
rhol = 0.1;pl = 0.1;ul = 15*sqrt(gamma*pl/rhol);
rhor = 0.1;pr = 0.1;ur =-15*sqrt(gamma*pr/rhor);

xpos = -.5;
for i = 1: idim
    xpos = xpos + 1.0/idim;
    
    if(xpos < 0)
        U(i,1) = rhol;
        U(i,2) = rhol*ul;
        U(i,3) = rhol*(0.5*ul^2+pl/(gamma-1)/rhol);
    else
        U(i,1) = rhor;
        U(i,2) = rhor*ur;
        U(i,3) = rhor*(0.5*ur^2+pr/(gamma-1)/rhor);
    end
end

end

% moving discontinuity
function U = test3(U,idim, gamma)
rhor = 10.0;pr = 1.0;ur =0.3*sqrt(gamma*pr/rhor);
rhol = 0.125;pl = 1.0;ul = 0.3*sqrt(gamma*pr/rhor);

xpos = -.5;
for i = 1: idim
    xpos = xpos + 1.0/idim;
    
    if(xpos < 0)
        U(i,1) = rhol;
        U(i,2) = rhol*ul;
        U(i,3) = rhol*(0.5*ul^2+pl/(gamma-1)/rhol);
    else
        U(i,1) = rhor;
        U(i,2) = rhor*ur;
        U(i,3) = rhor*(0.5*ur^2+pr/(gamma-1)/rhor);
    end
end

end

% strong expansion fan
function U = test4(U,idim, gamma)
%rhol = 1.0;pl = 2.0;ul = -2.5*sqrt(gamma*pl/rhol);
%rhor = 1.0;pr = 0.5;ur =  2.5*sqrt(gamma*pl/rhol);
rhol = 1.0;pl = 2.0;ul = -1000;
rhor = 1.0;pr = 2.0;ur =  1000;

xpos = -.5;
for i = 1: idim
    xpos = xpos + 1.0/idim;
    
    if(xpos < 0)
        U(i,1) = rhol;
        U(i,2) = rhol*ul;
        U(i,3) = rhol*(0.5*ul^2+pl/(gamma-1)/rhol);
    else
        U(i,1) = rhor;
        U(i,2) = rhor*ur;
        U(i,3) = rhor*(0.5*ur^2+pr/(gamma-1)/rhor);
    end
end

end

function U = BoundaryCorrection(U)
% Periodic
U(1,:)  = U(end-1,:);
U(end,:)= U(2,:);

% Periodic - 2 layers
%U(1,:)  = U(end-3,:);
%U(2,:)  = U(end-2,:);
%U(end,:)= U(4,:);
%U(end-1,:)= U(3,:);

% Reflective - X
%U(1,:) = U(2,:);
%U(end,:) = U(end-1,:);

%U(1,2) = -U(2,2);
%U(end,2) = -U(end-1,2);

% Transmissive
%{
rho(1)  = rho(2);
rho(end)= rho(end-1);

rhou(1)  = rhou(2);
rhou(end)= rhou(end-1);

e(1)  = e(2);
e(end)= e(end-1);
%}
end

function dt = cflcon(cflcoe, u, a, dx, t, tf)
dt = cflcoe * dx / max(abs(u) + a);
mustBeNonnegative(dt);

if((t+dt) > tf)
    dt = tf - t;
end

end

function U = update(U, flux, dt, dx)
U(2:end-1,:) =  U(2:end-1,:) - dt/dx*flux(2:end-1,:);

end

function flux = fluxes(flux,U,dx,dt,gamma,nvars)
rhol = U(1:end-1,1);
rhor = U(2:end,1);
ul = U(1:end-1,2)./rhol;
ur = U(2:end,2)  ./rhor;
pl = (U(1:end-1,3) - 0.5.*rhol.*ul.^2)*(gamma-1);
pr = (U(2:end,3)   - 0.5.*rhor.*ur.^2)*(gamma-1);

% Call riemann solver for inviscid burger
riemann_solver_t = @(rhol,ul,pl,rhor,ur,pr)(riemann_solver(rhol,ul,pl,rhor,ur,pr,gamma));
[~, ~, sampler] = arrayfun(riemann_solver_t,rhol,ul,pl,rhor,ur,pr, 'UniformOutput',false);

% Compute Godunov intercell flux
%flux = 0.5*ustar.*ustar;
for i = 1:numel(sampler)
    [rs,us,ps] = feval(sampler{i},0,dt);
    flux(i,:) = compute_flux([rs us ps],gamma);
end

end

function FU = compute_flux(W,gamma)
FU = [W(:,1).*W(:,2) W(:,1).*W(:,2).^2+W(:,3) W(:,2).*(W(:,3) + W(:,1).*(0.5*W(:,2).^2 + W(:,3)./W(:,1)/(gamma-1)))];
end

function ustar = riemann(ul,ur)
if(ul > ur)
    S = 0.5*(ul+ur);
    
    if(S >= 0.0)
        ustar = ul;
    else
        ustar = ur;
    end
else
    if(ul >= 0.0)
        ustar = ul;
    end
    
    if(ur <= 0.0)
        ustar = ur;
    end
    
    if(ul <= 0.0 && ur >= 0.0)
        ustar = 0.0;
    end
end

end

function phi = minmod3(a,b,c,t)
M1 = (a > 0) .* min(min(t*a,b),t*c);
M2 = (a < 0) .* max(max(t*a,b),t*c);
phi = (M1 + M2).*(a.*b.*c > 0);
end

function flux = Rusanov_fluxes(flux,U,dx,dt,gamma,nvars)

W = ToPrimitiveVars(U,gamma);
Wm1 = circshift(W, 1,1); %i-1
Wp1 = circshift(W,-1,1); %i+1


WLL = W + 0.5*(W - Wm1).*max(0,min(1,(Wp1-W)./(W - Wm1)));
WLL = circshift(WLL,1,1);
WLR = W + 0.5*(W - Wp1).*max(0,min(1,(Wm1-W)./(W - Wp1)));


ULL = ToConservedVars(WLL,gamma);
ULR = ToConservedVars(WLR,gamma);
ALL = (getSoundSpeed(ULL, gamma));
ALR = (getSoundSpeed(ULR, gamma));
a = max(abs(WLL(:,2))+ALL,abs(WLR(:,2))+ALR);
%flux = 0.5*(compute_flux(WLR,gamma)+compute_flux(WLL,gamma) - a.*(ULR-ULL));
%CM = (ALL+ALR)/2;
FP = euler_jacobian_flux_splitting(WLL,gamma, 1);
FN = euler_jacobian_flux_splitting(WLR,gamma,-1);
flux = FP + FN;
flux = flux(2:end,:);
flux = diff(flux);
flux = [zeros(1,3);flux;zeros(1,3)];
mustBeReal(flux);

end

function F = euler_jacobian_flux_splitting(W,gamma,sgn)
rho = W(:,1);
u = W(:,2);
c = getSoundSpeed(W, gamma);
L1 = (u + sgn*abs(u))/2;
L2 = (u+c + sgn*abs(u+c))/2;
L3 = (u-c + sgn*abs(u-c))/2;

F(:,1) = 2*(gamma-1).*L1 + L2 + L3;
F(:,2) = 2*(gamma-1).*L1.*u + L2.*(u + c) + L3.*(u - c);
F(:,3) = (gamma - 1) .* L1 .* u.^2 + L2.*(u + c).^2/2 + L3.*(u - c).^2/2 + (3-gamma)*(L2+L3).*c.^2/2/(gamma-1);

F = F .* rho/2/gamma;

end

function flux = AUSMDV_fluxes(flux,U,dx,dt,gamma,nvars)
flux = AUSM_SOLVER1D.AUSM_fluxes(flux,U,dx,dt,gamma,3);
flux = diff(flux);
flux = [zeros(1,3);flux;zeros(1,3)];
mustBeReal(flux);

end