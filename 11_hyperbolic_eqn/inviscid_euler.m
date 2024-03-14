clear; clc; close all;

%% Inviscid Burgers
gamma = 1.4;
t0 = 0;
tf = 0.2;
idim = 300;
ti = 0.002;
flux = zeros(idim+1,3);

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


while(t < tf)
    %U = BoundaryCorrection(U);
    dt = cflcon(0.9, U(:,2)./U(:,1), getSoundSpeed(U, gamma), dx, t, tif);
    t = t + dt;
    flux = fluxes(flux, U, dt, gamma);
    U = update(U, flux, dt, dx);
    
    if(t == tif)
        sol.rho = [sol.rho U(:,1)];
        sol.u = [sol.u U(:,2)./U(:,1)];
        sol.p = [sol.p getPressure(U,gamma)];
        sol.e = [sol.e getInternalEnergy(U)];
        nti = nti+1;
        tif = ti*nti;
        
        disp(t)
    end

end

%% Plot
field = 'rho';
surf(sol.(field).','EdgeColor','none');
view(2), axis tight;
colormap gray;

%% Plot 2

subplot(2,2,1)
plot(sol.('rho')(:,end),'.-');
subplot(2,2,2)
plot(sol.('u')(:,end),'.-');
subplot(2,2,3)
plot(sol.('p')(:,end),'.-');
subplot(2,2,4)
plot(sol.('e')(:,end),'.-');

%% Helper functions
function u = getVelU(U)
    u = U(:,2)./U(:,1);
end

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
    a = sqrt(gamma*getPressure(U,gamma)./U(:,1));
end

function U = test1(U,idim, gamma)
rhol = 1.0;ul = 0.75;pl = 1.0;
rhor = 0.125;ur = 0.0;pr = 0.1;

xpos = -0.5;
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
%{
rho(1)  = rho(end-1);
rho(end)= rho(2);

rhou(1)  = rhou(end-1);
rhou(end)= rhou(2);

e(1)  = e(end-1);
e(end)= e(2);
%}

% Reflective - X
U(1,:) = U(2,:);
U(end,:) = U(end-1,:);

U(1,2) = -U(2,2);
U(end,2) = -U(end-1,2);

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

if((t+dt) > tf)
    dt = tf - t;
end

end

function U = update(U, flux, dt, dx)
U(2:end-1,:) =  U(2:end-1,:) - dt/dx*diff(flux);

end

function flux = fluxes(flux, U, dt, gamma)
rhol = U(1:end-1,1);
rhor = U(2:end,1);
ul = U(1:end-1,2)./rhol;
ur = U(2:end,2)  ./rhor;
pl = (U(1:end-1,3) - 0.5.*rhol.*ul.^2)*(gamma-1);
pr = (U(2:end,3)   - 0.5.*rhor.*ur.^2)*(gamma-1);

%{
% Call riemann solver for inviscid burger
riemann_solver_t = @(rhol,ul,pl,rhor,ur,pr)(riemann_solver(rhol,ul,pl,rhor,ur,pr,gamma));
[~, ~, sampler] = arrayfun(riemann_solver_t,rhol,ul,pl,rhor,ur,pr, 'UniformOutput',false);

% Compute Godunov intercell flux
%flux = 0.5*ustar.*ustar;
for i = 1:numel(sampler)
    [rs,us,ps] = feval(sampler{i},0,dt);
    flux(i,:) = [rs*us, rs*us^2+ps, us*(ps + rs*(0.5*us^2 + ps/rs/(gamma-1)))];
end
%}

riemann_solver_t = @(rhol,ul,pl,rhor,ur,pr)(roe_solver(rhol,ul,pl,rhor,ur,pr,gamma));
[flux(:,1),flux(:,2),flux(:,3)] = arrayfun(riemann_solver_t,...
    rhol,ul,pl,rhor,ur,pr);

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