function [F1,F2,F3] = roe_solver(RHOL,UL,PL,RHOR,UR,PR,GAMMA)
%% [USTAR, PSTAR, SAMPLER] = riemann_solver(RHOL,UL,PL,RHOR,UR,PR,GAMMA)
%

HL = 0.5*UL^2 + GAMMA/(GAMMA-1)*PL/RHOL;
HR = 0.5*UR^2 + GAMMA/(GAMMA-1)*PR/RHOR;

% Compute the Roe average values
ubar = (sqrt(RHOL)*UL + sqrt(RHOR)*UR)/(sqrt(RHOL) + sqrt(RHOR));
hbar = (sqrt(RHOL)*HL + sqrt(RHOR)*HR)/(sqrt(RHOL) + sqrt(RHOR));
abar = sqrt((GAMMA-1)*max(0,hbar - 0.5*ubar^2));

% Compute the averaged eigenvalues lbar
lbar = [ubar - abar ubar ubar + abar];

% Compute the averaged right eigenvectors K
K = [1 1 1;ubar-abar ubar ubar+abar;hbar-ubar*abar 0.5*ubar^2 hbar+ubar*abar];

% Solve the wave strength coeffcients ai
du = [RHOR;RHOR*UR;0.5*RHOR*UR^2 + PR/(GAMMA-1)]-[RHOL;RHOL*UL;0.5*RHOL*UL^2 + PL/(GAMMA-1)];
ai = K\du;

% Compute F_i+1/2
flux = 0.5*(compute_flux([RHOL UL PL],GAMMA)+compute_flux([RHOR UR PR],GAMMA))...
    -0.5*(abs(lbar).*K)*ai;

% Entropy fix
%{
RSL = RHOL+ai(1);
USL = (RHOL*UL + ai(1)*(ubar-abar))/(RHOL+ai(1));
PSL = max(0,(GAMMA-1)*((0.5*RHOL*UL^2+PL/(GAMMA-1))+ai(1)*(hbar-ubar*abar)-0.5*RSL*USL^2));

RSR = RHOR-ai(end);
USR = (RHOR*UR - ai(end)*(ubar+abar))/(RHOR-ai(end));
PSR = max(0,(GAMMA-1)*((0.5*RHOR*UR^2+PR/(GAMMA-1))-ai(end)*(hbar+ubar*abar)-0.5*RSR*USR^2));

ASL = sqrt(GAMMA*PSL/RSL);
ASR = sqrt(GAMMA*PSR/RSR);
mustBeReal(ASL);
mustBeReal(ASR);

SL = UL - sqrt(GAMMA*PL/RHOL);
SR = USL - ASL;
if(SL < 0 && SR > 0) % Left Transonic Rarefaction
    lbar_fixed = SL*(SR - lbar(1))/(SR - SL);
    flux = compute_flux([RHOL UL PL],GAMMA) +lbar_fixed*ai(1).*K(:,1);
else
    SL = USR + ASR;
    SR = UR + sqrt(GAMMA*PR/RHOR);
    if(SL < 0 && SR > 0) % Right Transonic Rarefaction
        lbar_fixed = SR*(lbar(end)-SL)/(SR - SL);
        flux = compute_flux([RHOR UR PR],GAMMA) -lbar_fixed*ai(end).*K(:,end);
    end
end
%}


% Assign to output vars
F1 = flux(1);
F2 = flux(2);
F3 = flux(3);

end

%% HELPER FUNCTIONS
function FU = compute_flux(W,gamma)
FU = [W(:,1).*W(:,2); ... % rho*u
      W(:,1).*W(:,2).^2+W(:,3); ... % rho*u^2 + p
      W(:,2).*(W(:,3) + W(:,1).*(0.5*W(:,2).^2 + W(:,3)./W(:,1)/(gamma-1)))];
end

function pm = guessp(rhol,ul,pl,rhor,ur,pr,gamma)
quser = 2.0;
cl = sqrt(gamma*pl/rhol);
cr = sqrt(gamma*pr/rhor);

g_1 = (gamma-1)/2/gamma;
g_2 = 2*gamma/(gamma-1);
g_5 = 2.0/(gamma + 1.0);
g_6 = (gamma - 1.0)/(gamma + 1.0);

cup = 0.25*(rhol + rhor)*(cl + cr);
ppv = 0.5*(pl + pr) + 0.5*(ul - ur)*cup;
ppv = max(0,ppv);
pmin = min(pl, pr);
pmax = max(pl, pr);
qmax = pmax/pmin;

if(qmax <= quser && (pmin <= ppv && ppv <= pmax))
    % select PVRS Riemann solver
    pm = ppv;
else
    if(ppv < pmin)
        % Select Two-Rarefaction Riemann solver
        pm = ((cl+cr-0.5*(gamma-1)*(ur-ul))/(cl/pl^g_1 + cr/pr^g_1))^g_2;
    else
        % Select Two-Shock Riemann solver with PVRS as estimate
        gel = sqrt((g_5/rhol)/(g_6*pl + ppv));
        ger = sqrt((g_5/rhor)/(g_6*pr + ppv));
        pm = (gel*pl + ger*pr - (ur - ul))/(gel + ger);
    end
end

end