function [F1,F2,F3,F5] = roe_solver2d(RHOL,UL,VL,PL,RHOR,UR,VR,PR,GAMMA)
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

% Assign to output vars
%F1 = 0.5*(RHOL*UL+RHOR*UR)          -0.5*(abs(lm1)*K11*a1+abs(lm2)*K12*a2+abs(lm3)*K13*a3+abs(lm5)*K15*a5);
%F2 = 0.5*(RHOL*UL^2+PL+RHOR*UR^2+PR)-0.5*(abs(lm1)*K21*a1+abs(lm2)*K22*a2+abs(lm3)*K23*a3+abs(lm5)*K25*a5);
%F3 = 0.5*(RHOL*UL*VL+RHOR*UR*VR)    -0.5*(abs(lm1)*K31*a1+abs(lm2)*K32*a2+abs(lm3)*K33*a3+abs(lm5)*K35*a5);
%F5 = 0.5*(FXL+FXR)                  -0.5*(abs(lm1)*K51*a1+abs(lm2)*K52*a2+abs(lm3)*K53*a3+abs(lm5)*K55*a5);


% Entropy fix
RSL = RHOL+a1;
USL = (RHOL*UL + a1*(ubar-abar))/(RHOL+a1);
PSL = max(0,(GAMMA-1)*((RHOL*VLSQR+PL/(GAMMA-1))+a1*(hbar-ubar*abar)-0.5*RSL*USL^2));

RSR = RHOR-a5;
USR = (RHOR*UR - a5*(ubar+abar))/(RHOR-a5);
PSR = max(0,(GAMMA-1)*((RHOR*VRSQR+PR/(GAMMA-1))-a5*(hbar+ubar*abar)-0.5*RSR*USR^2));

ASL = sqrt(GAMMA*PSL/RSL);
ASR = sqrt(GAMMA*PSR/RSR);

SL = UL - sqrt(GAMMA*PL/RHOL);
SR = USL - ASL;
if(SL < 0 && SR > 0) % Left Transonic Rarefaction
        lm1fixed = SL*(SR - lm1)/(SR - SL);
        F1 = RHOL*UL + lm1fixed*K11*a1;
        F2 = RHOL*UL^2+PL + lm1fixed*K21*a1;
        F3 = RHOL*UL*VL + lm1fixed*K31*a1;
        F5 = FXL + lm1fixed*K51*a1;
else
    SL = USR + ASR;
    SR = UR + sqrt(GAMMA*PR/RHOR);
    if(SL < 0 && SR > 0) % Right Transonic Rarefaction
        lm5fixed = SR*(lm5-SL)/(SR - SL);
        F1 = RHOR*UR - lm5fixed*K15*a5;
        F2 = RHOR*UR^2+PR - lm5fixed*K25*a5;
        F3 = RHOR*UR*VR - lm5fixed*K35*a5;
        F5 = FXR - lm5fixed*K55*a5;
    
    else % Subject to normal solution
        % Assign to output vars
        F1 = 0.5*(RHOL*UL+RHOR*UR) -0.5*(abs(lm1)*K11*a1+abs(lm2)*K12*a2+abs(lm3)*K13*a3+abs(lm5)*K15*a5);
        F2 = 0.5*(RHOL*UL^2+PL+RHOR*UR^2+PR) -0.5*(abs(lm1)*K21*a1+abs(lm2)*K22*a2+abs(lm3)*K23*a3+abs(lm5)*K25*a5);
        F3 = 0.5*(RHOL*UL*VL+RHOR*UR*VR) -0.5*(abs(lm1)*K31*a1+abs(lm2)*K32*a2+abs(lm3)*K33*a3+abs(lm5)*K35*a5);
        F5 = 0.5*(FXL+FXR) -0.5*(abs(lm1)*K51*a1+abs(lm2)*K52*a2+abs(lm3)*K53*a3+abs(lm5)*K55*a5);

    end
end


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