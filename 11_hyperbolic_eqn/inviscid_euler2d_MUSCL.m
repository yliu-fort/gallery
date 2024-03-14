clear; clc; close all;
addpath Colormaps

%% Inviscid Burgers
gamma = 1.4;
t0 = 0;
tf = 2.0/10;
NX = 800;
NY = 400;
ti = 0.002;
gpuFlag = true;
if(gpuFlag)
    flux = gpuArray.zeros(NY+1,NX+1,4);
    FX = gpuArray.zeros(NY+1,NX+1,4);
    FY = gpuArray.zeros(NY+1,NX+1,4);
    U = gpuArray.zeros(NY+2,NX+2,4);
else
    flux = zeros(NY+1,NX+1,4);
    FX = zeros(NY+1,NX+1,4);
    FY = zeros(NY+1,NX+1,4);
    U = zeros(NY+2,NX+2,4);
end

dx = 0.005;
dy = 0.005;
[X,Y] = meshgrid([0:(NX+1)]*dy,[0:(NY+1)]*dx);

U = test3(U,dx,dy,NY+2,NX+2, gamma,NX/NY);

[ibcoeff,IA,IC,p1,p2,p3,p4,a1,a2,a3,a4,~] = ImmersedBoundary2D(X,Y,dx);

t = t0;
nti = 1;
tif = ti*nti;

sol.gamma = gamma;
sol.t0 = t0;
sol.NX = NX;
sol.NY = NY;
sol.ti = ti;
sol.gpuFlag = gpuFlag;
sol.U0 = gather(U);
sol.X = gather(X);
sol.Y = gather(Y);
sol.dx = dx;
sol.dy = dy;
sol.ib = {ibcoeff,IA,IC,p1,p2,p3,p4,a1,a2,a3,a4};

%% Main program
field = 'e';
h = mesh([],'EdgeColor','none','FaceColor','interp');
view(2),axis tight,colormap viridis%,caxis([-0.25 0]);
colorbar,axis equal
set(gca,'FontName','Times New Roman','FontSize',16);
set(gcf,'color','w','Position',[50 50 800 400])
drawnow
iter = 0;

while(t < tf)
    dt = cflcon(0.33, U, gamma, dx, dy, t, tif);
    t = t + dt;
    iter = iter + 1;
    
    % Advance
    U = BoundaryCorrection(U);
    U = IBBoundaryCorrection(U,ibcoeff,IA,IC,p1,p2,p3,p4,a1,a2,a3,a4,gamma,X,Y);
    try
        %[FX,FY] = KT_fluxes(FX,FY,U,dx,dy,dt,gamma);
        [FX,FY] = MUSCL_fluxes(FX,FY,U,dx,dy,dt,gamma);
        %[FX,FY] = MUSCL_fluxes_legacy(FX,FY,U,dx,dy,dt,gamma);
        %[FX,FY] = CIR_fluxes(FX,FY,U,dx,dy,dt,gamma);
        %[FX,FY] = AUSM_SOLVER2D.AUSMDV_fluxes(FX,FY,U,dx,dy,dt,gamma);
        %[FX,FY] = AUSM_SOLVER2D.MUSCL_AUSM_fluxes(FX,FY,U,dx,dy,dt,gamma);
        
        cflcon(1.0, update(U, FX, FY, dt, dx, dy), gamma, dx, dy, t, tif);
    catch ME
        disp("Solve failed, switch to 2nd stable solver.")

        [FX,FY] = KT_fluxes(FX,FY,U,dx,dy,dt,gamma);
        %[FX,FY] = MUSCL_fluxes_stable(FX,FY,U,dx,dy,dt,gamma);

    end
    
    U = update(U, FX, FY, dt, dx, dy);
    
    %     W = ToPrimitiveVars(U,gamma);
    %     for k = 1:size(U,3)
    %         tW = W(:,:,k);
    %         tW(IA) = 1e-20;
    %         tW(IC) = 1e-20;
    %         W(:,:,k)=tW;
    %     end
    %     U = ToConservedVars(W,gamma);
    
    if(t == tif)
        sol.t{nti} = t;
        tU = gather(U);
        sol.rho{nti} = tU(:,:,1);
        %sol.u{nti} = tU(:,:,2)./tU(:,:,1);
        %sol.v{nti} = tU(:,:,3)./tU(:,:,1);
        sol.m{nti} = sqrt(tU(:,:,2).^2 + tU(:,:,3).^2)./getSoundSpeed(tU,gamma);
        sol.p{nti} = getPressure(tU,gamma);
        sol.e{nti} = getInternalEnergy(tU);
        [drx,dry] = gradient(sol.rho{nti});
        sol.drho{nti} = max(abs(drx),abs(dry));
        
        %if(nti == 1),h.ZData = 0*sol.(field){1};end
        h.ZData = sol.(field){nti};
        drawnow
        
        nti = nti+1;
        tif = ti*nti;
        disp(t +": "+iter)
    end
    
end

sol.Uf = gather(U);
sol.tf = t;

%save('result.mat', 'sol','-v7.3');

%% Animate 1
field = 'rho';
h = mesh(X,Y,sol.(field){1},'EdgeColor','none','FaceColor','interp');
view(2),axis tight,colormap viridis,caxis([0 4]);
cb = colorbar('North','Color','w'),ylabel(cb,'Rho')
axis equal
set(gca,'FontName','Times New Roman','FontSize',12,'xticklabels',[],'yticklabels',[]);
set(gca, 'Position', [0,0,1,1])
set(gcf,'OuterPosition',[10 10 1366 768]);

for i = 1:numel(sol.(field))
    h.ZData = sol.(field){i};
    drawnow
end

%% Animate 2
field = 'm';
h = mesh(X,Y,sol.(field){1},'EdgeColor','none','FaceColor','flat');
view(2),axis tight,colormap coldhot,caxis([0 2]);
cb = colorbar('North','Color','k'),ylabel(cb,'Mach')
axis equal
set(gca,'FontName','Times New Roman','FontSize',12,'xticklabels',[],'yticklabels',[]);
set(gca, 'Position', [0,0,1,1])
set(gcf,'OuterPosition',[10 10 1366 768]);
for i = 1:numel(sol.(field))
    h.ZData = sol.(field){i};
    drawnow
end

%% Animate 3
field = 'e';
h = mesh(X,Y,sol.(field){1},'EdgeColor','none','FaceColor','flat');
view(2),axis tight,colormap coldhot%,caxis([0 2]);
cb = colorbar('North','Color','k'),ylabel(cb,'Temperature')
axis equal
set(gca,'FontName','Times New Roman','FontSize',12,'xticklabels',[],'yticklabels',[]);
set(gca, 'Position', [0,0,1,1])
set(gcf,'OuterPosition',[10 10 1366 768]);
for i = 1:numel(sol.(field))
    h.ZData = sol.(field){i};
    drawnow
end

%% Animate 4
field = 'drho';
h = mesh(X,Y,sol.(field){1},'EdgeColor','none','FaceColor','interp');
view(2),axis tight,colormap gray,caxis([-.1 0]);
cb = colorbar('North'),ylabel(cb,'$\nabla \rho$','interpreter','latex')
axis equal
set(gca,'FontName','Times New Roman','FontSize',12,'xticklabels',[],'yticklabels',[]);
set(gca, 'Position', [0,0,1,1])
set(gcf,'OuterPosition',[10 10 1366 768]);
for i = 1:numel(sol.(field))
    h.ZData = -sol.(field){i};
    drawnow
end


%% Init
% Blast
function U = test_carbuncle(U,dx,dy,nx,ny, gamma,aspect)
rhol = 1.0;ul = 0.0;vl = 0.0;pl = 2.5;
rhor = 0.125;ur = 0.0;vr = 0.0;pr = 0.25;

[xpos, ypos] = meshgrid(0.7 + linspace(-1,1+2*(aspect-1),ny),...
    0.7 + linspace(-1,1,nx));

T = xpos < 0;
%T = 1-(tanh(((xpos.^2 + ypos.^2) - 0.25^2)/(dx^2+dy^2)/64^2)+1)/2;
%T = (xpos.^2 + ypos.^2) < 0.2 & (xpos.^2 + ypos.^2) > 0.1;
%T = (xpos < -0.2) & (ypos < -0.2) & (xpos > -0.5) & (ypos > -0.5);
%T = xpos+0.001*ypos < 0;
U(:,:,1) = rhol.*T + rhor.*(1-T);
U(:,:,2) = rhol.*ul.*T + rhor.*ur.*(1-T);
U(:,:,3) = rhol.*vl.*T + rhor.*vr.*(1-T);
U(:,:,4) = (rhol.*0.5*(ul.^2 + vl.^2)+pl/(gamma-1)).*T...
    + (rhor.*0.5*(ur.^2 + vr.^2)+pr/(gamma-1)).*(1-T);

U(:,1:2:end,:) = U(:,1:2:end,:) + 1e-6;
U(:,2:2:end,:) = U(:,2:2:end,:) - 1e-6;
end

function U = test0(U,dx,dy,nx,ny, gamma,aspect)
rhol = 10;ul = 0.0;vl = 0.0;pl = 25;
rhor = 0.125;ur = 0.0;vr = 0.0;pr = 0.25;

[xpos, ypos] = meshgrid(0.7 + linspace(-1,1+2*(aspect-1),ny),...
    0.7 + linspace(-1,1,nx));

T = xpos > aspect;
%T = 1-(tanh(((xpos.^2 + ypos.^2) - 0.25^2)/(dx^2+dy^2)/64^2)+1)/2;
%T = (xpos.^2 + ypos.^2) < 0.2 & (xpos.^2 + ypos.^2) > 0.1;
%T = (xpos < -0.2) & (ypos < -0.2) & (xpos > -0.5) & (ypos > -0.5);
%T = xpos+0.001*ypos < 0;
U(:,:,1) = rhol.*T + rhor.*(1-T);
U(:,:,2) = rhol.*ul.*T + rhor.*ur.*(1-T);
U(:,:,3) = rhol.*vl.*T + rhor.*vr.*(1-T);
U(:,:,4) = (rhol.*0.5*(ul.^2 + vl.^2)+pl/(gamma-1)).*T...
    + (rhor.*0.5*(ur.^2 + vr.^2)+pr/(gamma-1)).*(1-T);

end

function U = test1(U,dx,dy,nx,ny, gamma,aspect)
rhol = 1.0;ul = 0.0;vl = 0.0;pl = 2.5;
rhor = 0.125;ur = 0.0;vr = 0.0;pr = 0.1;

[xpos, ypos] = meshgrid(0.7 + linspace(-1,1+2*(aspect-1),ny),...
    0.7 + linspace(-1,1,nx));

T = ((xpos-2.0).^2 + ypos.^2) < 0.15^2;
%T = 1-(tanh(((xpos.^2 + ypos.^2) - 0.25^2)/(dx^2+dy^2)/64^2)+1)/2;
%T = (xpos.^2 + ypos.^2) < 0.2 & (xpos.^2 + ypos.^2) > 0.1;
%T = (xpos < -0.2) & (ypos < -0.2) & (xpos > -0.5) & (ypos > -0.5);
%T = xpos+0.001*ypos < 0;
U(:,:,1) = rhol.*T + rhor.*(1-T);
U(:,:,2) = rhol.*ul.*T + rhor.*ur.*(1-T);
U(:,:,3) = rhol.*vl.*T + rhor.*vr.*(1-T);
U(:,:,4) = (rhol.*0.5*(ul.^2 + vl.^2)+pl/(gamma-1)).*T...
    + (rhor.*0.5*(ur.^2 + vr.^2)+pr/(gamma-1)).*(1-T);

end


function U = test2(U,dx,dy,nx,ny, gamma,aspect)
rhol = 1.0;ul = 0.0;vl = 0.0;pl = 100;
rhor = 1.0;ur = 0.0;vr = 0.0;pr = 1;

[xpos, ypos] = meshgrid(0.5 + linspace(-1,1+2*(aspect-1),ny),...
    0.7 + linspace(-1,1,nx));

%T = (xpos.^2 + ypos.^2) < 0.25^2;
T = xpos < 0;
U(:,:,1) = rhol.*T + rhor.*(1-T);
U(:,:,2) = rhol.*ul.*T + rhor.*ur.*(1-T);
U(:,:,3) = rhol.*vl.*T + rhor.*vr.*(1-T);
U(:,:,4) = (rhol.*0.5*(ul.^2 + vl.^2)+pl/(gamma-1)).*T...
    + (rhor.*0.5*(ur.^2 + vr.^2)+pr/(gamma-1)).*(1-T);


rho = 0.8;u = 0.0;v = 0.0;p = 1.0;

[xpos, ypos] = meshgrid(0.0 + linspace(-1,1+2*(aspect-1),ny),...
    0.07 + linspace(-1,1,nx));

T = (xpos.^2 + ypos.^2) < 0.25^2;

U(:,:,1) = rho.*T + U(:,:,1).*(1-T);
U(:,:,2) = rho.*u.*T + U(:,:,2).*(1-T);
U(:,:,3) = rho.*v.*T + U(:,:,3).*(1-T);
U(:,:,4) = (rho.*0.5*(u.^2 + v.^2)+p/(gamma-1)).*T...
    + U(:,:,4).*(1-T);

end

% Throttle
function U = test3(U,dx,dy,nx,ny, gamma,aspect)
rhol = 1.0;ul = 0.0;vl = 10.0;pl = 1;
rhor = 1.0;ur = 0.0;vr = 10.0;pr = 1;

[xpos, ypos] = meshgrid(0.7 + linspace(-1,1+2*(aspect-1),ny),...
    0.7 + linspace(-1,1,nx));

U(:,:,1) = rhor;
U(:,:,2) = rhor.*ur;
U(:,:,3) = rhor.*vr;
U(:,:,4) = (rhor.*0.5*(ur.^2 + vr.^2)+pr/(gamma-1));

U(:,1,1) = rhol;
U(:,1,2) = rhol.*ul;
U(:,1,3) = rhol.*vl;
U(:,1,4) = (rhol.*0.5*(ul.^2 + vl.^2)+pl/(gamma-1));

end
%% MUSCL_Hancock
function [FX,FY] = KT_fluxes(FX,FY,U,dx,dy,dt,gamma)
% CFL < 0.5 to eliminate strips
% Data reconstruction X
W = ToPrimitiveVars(U,gamma);
Wm1 = circshift(W, 1,1); %i-1
Wp1 = circshift(W,-1,1); %i+1

% WLL = W + 0.5*(W - Wm1).*max(0,min(1,(Wp1-W)./(W - Wm1)));
% WLL = circshift(WLL,1,1);
% WLR = W + 0.5*(W - Wp1).*max(0,min(1,(Wm1-W)./(W - Wp1)));

dW = mc(Wp1 - W, W - Wm1);

WLL = W + 0.5*dW;
WLL = circshift(WLL,1,1);
WLR = W - 0.5*dW;

ULL = ToConservedVars(WLL,gamma);
ULR = ToConservedVars(WLR,gamma);

ALL = getSoundSpeed(ULL, gamma);
ALR = getSoundSpeed(ULR, gamma);
a = max(abs(WLL(:,:,2))+ALL,abs(WLR(:,:,2))+ALR);
FXstar = 0.5*(compute_fluxX(WLR,gamma)+compute_fluxX(WLL,gamma) - a.*(ULR-ULL));


% Data reconstruction Y
Wm1 = circshift(W, 1,2); %i-1
Wp1 = circshift(W,-1,2); %i+1


% WLL = W + 0.5*(W - Wm1).*max(0,min(1,(Wp1-W)./(W - Wm1)));
% WLL = circshift(WLL,1,2);
% WLR = W + 0.5*(W - Wp1).*max(0,min(1,(Wm1-W)./(W - Wp1)));

dW = mc(Wp1 - W, W - Wm1);

WLL = W + 0.5*dW;
WLL = circshift(WLL,1,2);
WLR = W - 0.5*dW;

ULL = ToConservedVars(WLL,gamma);
ULR = ToConservedVars(WLR,gamma);


ALL = getSoundSpeed(ULL, gamma);
ALR = getSoundSpeed(ULR, gamma);
a = max(abs(WLL(:,:,3))+ALL,abs(WLR(:,:,3))+ALR);
FYstar = 0.5*(compute_fluxY(WLR,gamma)+compute_fluxY(WLL,gamma) - a.*(ULR-ULL));


FX(:,1:end-1,:) = FXstar(2:end,2:end-1,:);
FY(1:end-1,:,:) = FYstar(2:end-1,2:end,:);

end

function phi = minmod(a,b)
M1 = (a > 0) .* min(min((a+b)/2,2*a),2*b);
M2 = (a < 0) .* max(max((a+b)/2,2*a),2*b);
phi = (M1 + M2).*(a.*b > 0);
end

function phi = mc(a,b)
d = (a > 0) & (b > 0) -(a < 0) & (b < 0) ;
phi = d.*min(min(2*abs(a),2*abs(b)),abs(a+b)/2);
end

function [FX,FY] = MUSCL_fluxes(FX,FY,U,dx,dy,dt,gamma)
% Data reconstruction
W = ToPrimitiveVars(U,gamma);
dWx = circshift(W, -1,1) - W; % i+1 - i
dWy = circshift(W, -1,2) - W; % j+1 - j

dWxL = circshift(dWx,1,1); % i - i-1
dWxR = dWx; % i+1 - i
dWyL = circshift(dWy,1,2);
dWyR = dWy;

dWxr = minmod(dWxL,dWxR);
dWyr = minmod(dWyL,dWyR);
%dWxr = mc(dWxL,dWxR);
%dWyr = mc(dWyL,dWyR);

% PP limiter
Vmax =  1e-20;
Vmin = -1e-20;
for i = -1:1
    for j = -1:1
        Vmax = max(Vmax, circshift(circshift(W,i,1),j,2)-W);
        Vmin = min(Vmin, circshift(circshift(W,i,1),j,2)-W);
    end
end

V = 2*min(abs(Vmin),abs(Vmax))./( abs(dWxr) + abs(dWyr) );
dWxr = min(1,V).*dWxr;
dWyr = min(1,V).*dWyr;

% Evolution
WL = W - 0.5.*dWxr;WR = W + 0.5.*dWxr;
WD = W - 0.5.*dWyr;WU = W + 0.5.*dWyr;

fluxL = compute_fluxX(WL,gamma);
fluxR = compute_fluxX(WR,gamma);
fluxD = compute_fluxY(WD,gamma);
fluxU = compute_fluxY(WU,gamma);

% prediction
Ustar = U - 0.5*dt/dx.*(fluxR - fluxL) - 0.5*dt/dy.*(fluxU - fluxD);
Wstar = ToPrimitiveVars(Ustar,gamma);
WL = Wstar - 0.5.*dWxr;WR = Wstar + 0.5.*dWxr;
WD = Wstar - 0.5.*dWyr;WU = Wstar + 0.5.*dWyr;

[f1,f2,f3,f4] = arrayfun(@rhll_solver2d,...
    WR(1:end-1,2:end-1,1),WR(1:end-1,2:end-1,2),WR(1:end-1,2:end-1,3),WR(1:end-1,2:end-1,4),...
    WL(2:end  ,2:end-1,1),WL(2:end  ,2:end-1,2),WL(2:end  ,2:end-1,3),WL(2:end  ,2:end-1,4),...
    gamma+zeros(size(WR(1:end-1,2:end-1,1))));

FXstar = cat(3,f1,f2,f3,f4);

[f1,f3,f2,f4] = arrayfun(@rhll_solver2d,...
    WU(2:end-1,1:end-1,1),WU(2:end-1,1:end-1,3),WU(2:end-1,1:end-1,2),WU(2:end-1,1:end-1,4),...
    WD(2:end-1,2:end  ,1),WD(2:end-1,2:end  ,3),WD(2:end-1,2:end  ,2),WD(2:end-1,2:end  ,4),...
    gamma+zeros(size(WU(2:end-1,1:end-1,1))));

FYstar = cat(3,f1,f2,f3,f4);

FX(:,1:end-1,:) = FXstar;
FY(1:end-1,:,:) = FYstar;

%
end

function [FX,FY] = MUSCL_fluxes_legacy(FX,FY,U,dx,dy,dt,gamma)
% Data reconstruction
dUx = diff(U, 1, 1);
dUy = diff(U, 1, 2);

Zx = 0.*dUx(1,:,:);
Zy = 0.*dUy(:,1,:);

dl = [Zx;dUx(1:end-1,:,:);Zx];
dr = [Zx;dUx(2:end  ,:,:);Zx];
dd = [Zy dUy(:,1:end-1,:) Zy];
du = [Zy dUy(:,2:end  ,:) Zy];

% Minmod limiter
beta = 1;
dxi = (dr > 0).*max(0, max(min(beta*dl,dr), min(dl,beta*dr))) ...
    + (dr < 0).*min(0, min(max(beta*dl,dr), max(dl,beta*dr)));
dyi = (du > 0).*max(0, max(min(beta*dd,du), min(dd,beta*du))) ...
    + (du < 0).*min(0, min(max(beta*dd,du), max(dd,beta*du)));

% Evolution
UL = U - 0.5.*dxi;UR = U + 0.5.*dxi;
UD = U - 0.5.*dyi;UU = U + 0.5.*dyi;

fluxL = compute_fluxX(ToPrimitiveVars(UL, gamma),gamma);
fluxR = compute_fluxX(ToPrimitiveVars(UR, gamma),gamma);
fluxD = compute_fluxY(ToPrimitiveVars(UD, gamma),gamma);
fluxU = compute_fluxY(ToPrimitiveVars(UU, gamma),gamma);

predFlux = 0.5*dt/dx.*(fluxL - fluxR) + 0.5*dt/dy.*(fluxD - fluxU);
WL = ToPrimitiveVars(UL + predFlux,gamma);
WR = ToPrimitiveVars(UR + predFlux,gamma);
WD = ToPrimitiveVars(UD + predFlux,gamma);
WU = ToPrimitiveVars(UU + predFlux,gamma);

[f1,f2,f3,f4] = arrayfun(@rhll_solver2d,...
    WR(1:end-1,2:end-1,1),WR(1:end-1,2:end-1,2),WR(1:end-1,2:end-1,3),WR(1:end-1,2:end-1,4),...
    WL(2:end  ,2:end-1,1),WL(2:end  ,2:end-1,2),WL(2:end  ,2:end-1,3),WL(2:end  ,2:end-1,4),...
    gamma+zeros(size(WR(1:end-1,2:end-1,1))));

FX(:,1:end-1,:) = cat(3,f1,f2,f3,f4);

[f1,f3,f2,f4] = arrayfun(@rhll_solver2d,...
    WU(2:end-1,1:end-1,1),WU(2:end-1,1:end-1,3),WU(2:end-1,1:end-1,2),WU(2:end-1,1:end-1,4),...
    WD(2:end-1,2:end  ,1),WD(2:end-1,2:end  ,3),WD(2:end-1,2:end  ,2),WD(2:end-1,2:end  ,4),...
    gamma+zeros(size(WU(2:end-1,1:end-1,1))));

FY(1:end-1,:,:) = cat(3,f1,f2,f3,f4);

end

%% Boundary Correction

function U = IBBoundaryCorrection(U,ibcoeff,ind,oind,p1,p2,p3,p4,a1,a2,a3,a4,gamma,X,Y)
W = ToPrimitiveVars(U,gamma);
for k = 1:size(U,3)
    tW = W(:,:,k);
    tW(oind) = 0.1;
    Wip = (a1.*tW(p1)+a2.*tW(p2)+a3.*tW(p3)+a4.*tW(p4));
    if(k == 1 || k == size(U,3)) % Density | Pressure
        tW(ind) = Wip;
    else  % Momentum
        tW(ind) = (1-ibcoeff).*Wip;
        tW(oind) = 0;
    end
    
    W(:,:,k)=tW;
end

U = ToConservedVars(W,gamma);
end

function U = IBBoundaryCorrectionAlt(U,ibcoeff,ind,p1,p2,p3,p4,a1,a2,a3,a4,gamma)

for k = 1:size(U,3)
    tU = U(:,:,k);
    Uip = (a1.*tU(p1)+a2.*tU(p2)+a3.*tU(p3)+a4.*tU(p4));
    if(k == 1) % Density
        tU(ind) = Uip;
    else % Stagnation Pressure
        if (k == size(U,3))
            tU(ind) = Uip;
        else % momentum
            tU(ind) = (1-ibcoeff).*Uip;
        end
    end
    
    U(:,:,k)=tU;
end

end

function U = BoundaryCorrection(U)
% box reflection test
% U(200,200,:) = U(200,199,:);
% U(200,200,2) = -U(200,199,2);
% U(200,200,3) = -U(200,199,3);

% Periodic
% U(1,:,:) = U(end-1,:,:);
% U(end,:,:) = U(2,:,:);

%U(:,1,:) = U(:,end-1,:);
%U(:,end,:) = U(:,2,:);


% Reflective noslip
% U(1,:,:) = U(2,:,:);
% U(end,:,:) = U(end-1,:,:);
%
% U(1,:,2) = -U(2,:,2);
% U(end,:,2) = -U(end-1,:,2);
% U(1,:,3) = -U(2,:,3);
% U(end,:,3) = -U(end-1,:,3);


% U(:,1,:) = U(:,2,:);
% U(:,end,:) = U(:,end-1,:);
%
% U(:,1,2) = -U(:,2,2);
% U(:,end,2) = -U(:,end-1,2);
% U(:,1,3) = -U(:,2,3);
% U(:,end,3) = -U(:,end-1,3);


% Transmissive
U(1,:,:) = U(2,:,:);
U(end,:,:) = U(end-1,:,:);

%U(:,1,:) = U(:,2,:);
U(:,end,:) = U(:,end-1,:);

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

function [a,marker] = getSoundSpeedSafe(U, gamma)
por = getPressure(U,gamma)./U(:,:,1);
marker = por > 0;
por(~marker) = 0;
a = sqrt(gamma*por);
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

function dt = cflcon(cflcoe, U, gamma, dx, dy, t, tf)
a = getSoundSpeed(U,gamma);
S = reshape(min(dx./(abs(getVelU(U)) + a), dy./(abs(getVelV(U)) + a)),[],1);
dt = 0.5*cflcoe*( min( S ) );

if((t+dt) > tf)
    dt = tf - t;
end

dt = gather(dt);

end

function U = update(U, fluxX, fluxY, dt, dx, dy)
U(2:end-1,2:end-1,:) =  U(2:end-1,2:end-1,:) ...
    - dt/dx*diff(fluxX(:,1:end-1,:), 1, 1) ...
    - dt/dy*diff(fluxY(1:end-1,:,:), 1, 2);

end

function [FX,FY] = CIR_fluxes(FX,FY,U,dx,dy,dt,gamma)
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
FX(:,1:end-1,:) = compute_fluxX(cat(3,rs,us,vs,ps), gamma);

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
FY(1:end-1,:,:) = compute_fluxY(cat(3,rs,us,vs,ps), gamma);
end

function FU = compute_fluxX(W,gamma)
FU = cat(3,  W(:,:,1).*W(:,:,2) ...
    ,W(:,:,1).*W(:,:,2).^2+W(:,:,4)...
    ,W(:,:,1).*W(:,:,2).*W(:,:,3)...
    ,W(:,:,2).*(W(:,:,4) + W(:,:,1).*(0.5*(W(:,:,2).^2+W(:,:,3).^2) + W(:,:,4)./W(:,:,1)/(gamma-1))) );

end

function FV = compute_fluxY(W,gamma)
FV = cat(3,  W(:,:,1).*W(:,:,3) ...
    ,W(:,:,1).*W(:,:,2).*W(:,:,3)...
    ,W(:,:,1).*W(:,:,3).^2+W(:,:,4)...
    ,W(:,:,3).*(W(:,:,4) + W(:,:,1).*(0.5*(W(:,:,2).^2+W(:,:,3).^2) + W(:,:,4)./W(:,:,1)/(gamma-1))) );
end