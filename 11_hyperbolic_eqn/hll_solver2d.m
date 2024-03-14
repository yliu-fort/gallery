function [F1,F2,F3,F5] = hll_solver2d(RHOL,UL,VL,PL,RHOR,UR,VR,PR,GAMMA)
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
AL = sqrt(GAMMA*PL/RHOL);
AR = sqrt(GAMMA*PR/RHOR);
SL = min(0,min(UL-AL,ubar - abar));
SR = max(0,max(UR+AR,ubar + abar));

% Assign to output vars
F1 = compute_hllflux(SL,SR,RHOL*UL,RHOR*UR,RHOL,RHOR);
F2 = compute_hllflux(SL,SR,RHOL*UL^2+PL,RHOR*UR^2+PR,RHOL*UL,RHOR*UR);
F3 = compute_hllflux(SL,SR,RHOL*UL*VL,RHOR*UR*VR,RHOL*VL,RHOR*VR);
F5 = compute_hllflux(SL,SR,FXL,FXR,RHOL*VLSQR + PL/(GAMMA-1),RHOR*VRSQR + PR/(GAMMA-1));


end

function FHLL = compute_hllflux(SL,SR,FL,FR,UL,UR)
FHLL = (SR*FL - SL*FR + SL*SR*(UR-UL))/(SR-SL);
end