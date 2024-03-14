clear; clc; close all;

%% Inviscid Burgers
gamma = 1.4;
t0 = 0;
tf = 1.0;
idim = 512;
ti = 0.005;
flux = gpuArray.zeros(idim+1,idim+1,4);
FX = gpuArray.zeros(idim+1,idim+1,4);
FY = gpuArray.zeros(idim+1,idim+1,4);
U = gpuArray.zeros(idim+2,idim+2,4);

dx = 1.0/idim;
dy = 1.0/idim;

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
    dt = cflcon(0.9, U, gamma, dx, dy, t, tif);
    t = t + dt;
    
    % Scan X
    U = BoundaryCorrection(U,1);
    U = BoundaryCorrection(U,2);
    FX = fluxesX(FX, U, dt, gamma);
    FY = fluxesY(FY, U, dt, gamma);
    U = update(U, FX, FY, dt, dx, dy);
    
    % Scan Y
    %U = BoundaryCorrection(U,2);
    %flux = fluxesY(flux, U, 0.5*dt, gamma);
    %U = updateY(U, flux, 0.50*dt, dy);
    
    % Scan X
    %U = BoundaryCorrection(U,1);
    %flux = fluxesX(flux, U, 0.25*dt, gamma);
    %U = updateX(U, flux, 0.25*dt, dx);
    
    if(t == tif)
        tU = gather(U);
        sol.rho{nti} = tU(:,:,1);
        sol.u{nti} = tU(:,:,2)./tU(:,:,1);
        sol.v{nti} = tU(:,:,3)./tU(:,:,1);
        sol.p{nti} = getPressure(tU,gamma);
        sol.e{nti} = getInternalEnergy(tU);
        nti = nti+1;
        tif = ti*nti;
        
        disp(t)
    end
    
end

%% Animate
field = 'rho';
h = surf(sol.(field){1},'EdgeColor','none','FaceColor','interp');
view(2),axis tight,colormap jet,caxis([0 1]);

for i = 1:numel(sol.(field))
    h.ZData = sol.(field){i};
    drawnow
end


%% Helper functions
function tU = getVelU(U)
tU = U(:,:,2)./U(:,:,1);
end

function tV = getVelV(U)
tV = U(:,:,3)./U(:,:,1);
end

function umag = getVelMagSqr(U)
umag = sum((U(:,:,2:(end-1))./U(:,:,1)).^2,3);
end

function e = getInternalEnergy(U)
e = (U(:,:,end) - 0.5*U(:,:,1).*getVelMagSqr(U))./U(:,:,1);
end

function p = getPressure(U, gamma)
p = (gamma-1)*(U(:,:,end) - 0.5*U(:,:,1).*getVelMagSqr(U));
end

function a = getSoundSpeed(U, gamma)
a = sqrt(gamma*getPressure(U,gamma)./U(:,:,1));
end

function W = ToPrimitiveVars(U, gamma)
W(:,:,1) = U(:,:,1);
W(:,:,2) = U(:,:,2)./U(:,:,1);
W(:,:,3) = U(:,:,3)./U(:,:,1);
W(:,:,4) = (U(:,:,4) - 0.5.*(U(:,:,2).^2+U(:,:,3).^2)./U(:,:,1))*(gamma-1);
end

function U = ToConservedVars(W, gamma)
U(:,:,1) = W(:,:,1);
U(:,:,2) = W(:,:,1).*W(:,:,2);
U(:,:,3) = W(:,:,1).*W(:,:,3);
U(:,:,4) = 0.5*W(:,:,1).*(W(:,:,2).^2+W(:,:,3).^2) + W(:,:,4)/(gamma-1);
end

function U = test1(U,idim, gamma)
rhol = 1.0;ul = 0.0;vl = 0.0;pl = 1.0;
rhor = 0.125;ur = 0.0;vr = 0.0;pr = 0.1;

[xpos, ypos] = meshgrid(-0.5 + linspace(-1,1,idim),...
    -0.3 + linspace(-1,1,idim));

T = (xpos.^2 + ypos.^2) < 0.1;
U(:,:,1) = rhol.*T + rhor.*(1-T);
U(:,:,2) = rhol.*ul.*T + rhor.*ur.*(1-T);
U(:,:,2) = rhol.*vl.*T + rhor.*vr.*(1-T);
U(:,:,4) = (rhol.*(0.5*(ul.^2 + vl.^2)+pl/(gamma-1)./rhol)).*T...
         + (rhor.*(0.5*(ur.^2 + vr.^2)+pr/(gamma-1)./rhor)).*(1-T);

end

function U = BoundaryCorrection(U, dim)
% Periodic
%{
rho(1)  = rho(end-1);
rho(end)= rho(2);

rhou(1)  = rhou(end-1);
rhou(end)= rhou(2);

e(1)  = e(end-1);
e(end)= e(2);
%}

if(dim == 1) % Reflective slip - X
U(1,:,:) = U(2,:,:);
U(end,:,:) = U(end-1,:,:);

U(1,:,2) = -U(2,:,2);
U(end,:,2) = -U(end-1,:,2);
end

if(dim == 2) % Reflective slip - Y
U(:,1,:) = U(:,2,:);
U(:,end,:) = U(:,end-1,:);

U(:,1,3) = -U(:,2,3);
U(:,end,3) = -U(:,end-1,3);
end

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

function dt = cflcon(cflcoe, U, gamma, dx, dy, t, tf)
a = getSoundSpeed(U,gamma);
S = reshape((abs(getVelU(U)) + a)*dy + (abs(getVelV(U)) + a)*dx,[],1);
dt = cflcoe*( dx*dy / max( S ) );

if((t+dt) > tf)
    dt = tf - t;
end

end

function U = update(U, fluxX, fluxY, dt, dx, dy)
U(2:end-1,2:end-1,:) =  U(2:end-1,2:end-1,:) ...
    - dt/dx*diff(fluxX(:,1:end-1,:), 1, 1) ...
    - dt/dy*diff(fluxY(1:end-1,:,:), 1, 2);

end

function flux = fluxesX(flux, U, dt, gamma)
W = ToPrimitiveVars(U,gamma);
rhol = W(1:end-1,2:end-1,1);
rhor = W(2:end  ,2:end-1,1);
ul   = W(1:end-1,2:end-1,2);
ur   = W(2:end  ,2:end-1,2);
vl   = W(1:end-1,2:end-1,3);
vr   = W(2:end  ,2:end-1,3);
pl   = W(1:end-1,2:end-1,4);
pr   = W(2:end  ,2:end-1,4);

% Call riemann solver for inviscid burger
%riemann_solver_t = @(rhol,ul,vl,pl,rhor,ur,vr,pr)(riemann_solver_md(0,dt,rhol,ul,vl,pl,rhor,ur,vr,pr,gamma));
[rs,us,vs,ps] = arrayfun(@riemann_solver_md, ...
    0,dt,rhol,ul,vl,pl,rhor,ur,vr,pr,gamma);

% Compute Godunov intercell flux
flux(:,1:end-1,1) = rs.*us;
flux(:,1:end-1,2) = rs.*us.^2+ps;
flux(:,1:end-1,3) = rs.*us.*vs;
flux(:,1:end-1,4) = us.*(ps + rs.*(0.5*(us.^2+vs.^2) + ps./rs/(gamma-1)));

end

function flux = fluxesY(flux, U, dt, gamma)
W = ToPrimitiveVars(U,gamma);
rhol = W(2:end-1,1:end-1,1);
rhor = W(2:end-1,2:end,1);
ul   = W(2:end-1,1:end-1,2);
ur   = W(2:end-1,2:end,2);
vl   = W(2:end-1,1:end-1,3);
vr   = W(2:end-1,2:end,3);
pl   = W(2:end-1,1:end-1,4);
pr   = W(2:end-1,2:end,4);

% Call riemann solver for inviscid burger
%riemann_solver_t = @(rhol,ul,vl,pl,rhor,ur,vr,pr)(riemann_solver_md(0,dt,rhol,ul,vl,pl,rhor,ur,vr,pr,gamma));
[rs,vs,us,ps] = arrayfun(@riemann_solver_md, ...
    0,dt,rhol,vl,ul,pl,rhor,vr,ur,pr,gamma);

% Compute Godunov intercell flux
flux(1:end-1,:,1) = rs.*vs;
flux(1:end-1,:,2) = rs.*vs.*us;
flux(1:end-1,:,3) = rs.*vs.^2+ps;
flux(1:end-1,:,4) = vs.*(ps + rs.*(0.5*(us.^2+vs.^2) + ps./rs/(gamma-1)));

end