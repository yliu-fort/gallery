function [F1,F2,F3,F5] = rhll_solver2d(RHOL,UL,VL,PL,RHOR,UR,VR,PR,GAMMA)
%% [USTAR, PSTAR, SAMPLER] = riemann_solver(RHOL,UL,PL,RHOR,UR,PR,GAMMA)
%
VLSQR = 0.5*(UL^2+VL^2);
VRSQR = 0.5*(UR^2+VR^2);

FXL = UL*(RHOL*VLSQR + GAMMA*PL/(GAMMA-1));
FXR = UR*(RHOR*VRSQR + GAMMA*PR/(GAMMA-1));

HLB = (RHOL*VLSQR + GAMMA/(GAMMA-1)*PL)/sqrt(RHOL);
HRB = (RHOR*VRSQR + GAMMA/(GAMMA-1)*PR)/sqrt(RHOR);

% Compute the Roe average values
RHOA = sqrt(RHOL) + sqrt(RHOR);
ubar = (sqrt(RHOL)*UL + sqrt(RHOR)*UR)/RHOA;
vbar = (sqrt(RHOL)*VL + sqrt(RHOR)*VR)/RHOA;
hbar = (HLB + HRB)/RHOA;
abar = sqrt( (GAMMA-1)*(hbar - 0.5*(ubar^2+vbar^2)) );

% Compute the averaged eigenvalues lbar
lm1 = ubar - abar;
lm2 = ubar;
lm3 = ubar;
lm5 = ubar + abar;

% Compute the averaged right eigenvectors K
K11 = 1;K12 = 1;K13 = 0;K15 = 1;
K21 = ubar-abar;K22 = ubar;K23 = 0;K25 = ubar+abar;
K31 = vbar; K32 = vbar; K33 = 1; K35 = vbar;
K51 = hbar-ubar*abar; K52 = 0.5*(ubar^2+vbar^2); K53 = vbar; K55 = hbar+ubar*abar;

du1 = RHOR-RHOL;
du2 = RHOR*UR-RHOL*UL;
du3 = RHOR*VR-RHOL*VL;
du5 = (RHOR*VRSQR-RHOL*VLSQR) + (PR-PL)/(GAMMA-1);
du5 = du5 - (du3 - vbar*du1)*vbar;

% Eave strength coeffcients ai
a2 = (GAMMA-1)/abar^2*(du1*(hbar-ubar^2) + ubar*du2 - du5);
a1 = 1/2/abar*(du1*(ubar+abar) - du2 - abar*a2);
a3 = du3 - vbar*du1;
a5 = du1 - (a1 + a2);

% HLL::Compute the fastest signal velocities
AL = sqrt(GAMMA*PL/RHOL);
AR = sqrt(GAMMA*PR/RHOR);
SL = min(0,min(UL-AL,ubar - abar));
SR = max(0,max(UR+AR,ubar + abar));

% Rotate::Compute co1 and co2
qn = sqrt((UR-UL)^2+(VR-VL)^2);
co1 = 0;
co2 = 1;
if(qn > 1e-6)
    co1 = abs(UR-UL)/qn; % HLL
    co2 = abs(VR-VL)/qn; % Roe
end

% RHLL::Modify eigenvalues lbar
C1 = co2*(SR+SL)/(SR-SL);
C2 = 2*co1*SR*SL/(SR-SL);

lm1f = co2*abs(lm1) - (C1*lm1 + C2);
lm2f = co2*abs(lm2) - (C1*lm2 + C2);
lm3f = co2*abs(lm3) - (C1*lm3 + C2);
lm5f = co2*abs(lm5) - (C1*lm5 + C2);

% Assign to output vars
%F1 = 0.5*(RHOL*UL+RHOR*UR) -0.5*((lm1f)*K11*a1+(lm2f)*K12*a2+(lm3f)*K13*a3+(lm5f)*K15*a5);
%F2 = 0.5*(RHOL*UL^2+PL+RHOR*UR^2+PR) -0.5*((lm1f)*K21*a1+(lm2f)*K22*a2+(lm3f)*K23*a3+(lm5f)*K25*a5);
%F3 = 0.5*(RHOL*UL*VL+RHOR*UR*VR) -0.5*((lm1f)*K31*a1+(lm2f)*K32*a2+(lm3f)*K33*a3+(lm5f)*K35*a5);
%F5 = 0.5*(FXL+FXR) -0.5*((lm1f)*K51*a1+(lm2f)*K52*a2+(lm3f)*K53*a3+(lm5f)*K55*a5);

F1 = (SR*RHOL*UL-SL*RHOR*UR)/(SR-SL) -0.5*((lm1f)*K11*a1+(lm2f)*K12*a2+(lm3f)*K13*a3+(lm5f)*K15*a5);
F2 = (SR*(RHOL*UL^2+PL)-SL*(RHOR*UR^2+PR))/(SR-SL) -0.5*((lm1f)*K21*a1+(lm2f)*K22*a2+(lm3f)*K23*a3+(lm5f)*K25*a5);
F3 = (SR*RHOL*UL*VL-SL*RHOR*UR*VR)/(SR-SL) -0.5*((lm1f)*K31*a1+(lm2f)*K32*a2+(lm3f)*K33*a3+(lm5f)*K35*a5);
F5 = (SR*FXL-SL*FXR)/(SR-SL) -0.5*((lm1f)*K51*a1+(lm2f)*K52*a2+(lm3f)*K53*a3+(lm5f)*K55*a5);


end

%% HELPER FUNCTIONS
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

function FHLL = compute_hllflux(SL,SR,FL,FR,UL,UR)
FHLL = SR*SL/(SR-SL)*(UR-UL)-0.5*(SR+SL)/(SR-SL)*(FR-FL);
end