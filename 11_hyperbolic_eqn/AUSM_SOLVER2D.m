classdef AUSM_SOLVER2D
    methods(Static)
        function upm = u_splitting(u,c,sgn)
            upm = 0.5*(u+sgn.*abs(u));
            usub = sgn.*0.25./c.*(u+sgn.*c).^2;
            upm(abs(u)<=c) = usub(abs(u)<=c);
            %upm(~isfinite(upm))=0;
        end
        
        function ppm = p_splitting(p,u,upm,c,sgn)
            ppm = 1./u;
            psub = 1./c.*(sgn.*2-u./c);
            ppm(abs(u)<=c) = psub(abs(u)<=c);
            ppm = ppm .* p .* upm;
            %ppm(~isfinite(ppm))=0;
        end
        
        function upm = u_splittingd(u,cm,a,sgn)
            upm = 0.5*(u+sgn.*abs(u));
            usub = a.*(sgn*(u+sgn*cm).^2./4./cm) + (1-a).*(u+sgn*abs(u))/2;
            upm(abs(u)<=cm) = usub(abs(u)<=cm);
            %upm(~isfinite(upm))=0;
        end
        
        function ppm = p_splittingd(p,u,cm,sgn)
            ppm = p.*(1+sgn*sign(u))/2;
            psub = p.*(u./cm+sgn).^2/4.*(2-sgn*u./cm);
            ppm(abs(u)<=cm) = psub(abs(u)<=cm);
            %ppm(~isfinite(ppm))=0;
        end
        
        function marker = IsShockFixLoc(UR,VR,AR)
            UL = circshift(UR,1,1);
            AL = circshift(AR,1,1);
            markerX = ((UL-AL > 0) & (UR-AR < 0)) | ((UL+AL > 0) & (UR+AR < 0));
            
            VL = circshift(VR,1,2);
            AL = circshift(AR,1,2);
            markerY = ((VL-AL > 0) & (VR-AR < 0)) | ((VL+AL > 0) & (VR+AR < 0));
            
            marker = markerX | circshift(markerX,1,1) | markerY | circshift(markerY,1,2);
        end
        
        function [marker, XCOR,YCOR] = IsSonicPoint(UR,VR,AR)
            UL = circshift(UR,1,1);
            AL = circshift(AR,1,1);
            markerA = ((UL-AL < 0) & (UR-AR > 0));
            markerB = ((UL+AL < 0) & (UR+AR > 0));
            RCOR = markerA.*(UR-AR) + markerB.*(UR+AR);
            LCOR = markerA.*(UL-AL) + markerB.*(UL+AL);
            XCOR = RCOR-LCOR;
            XCOR(markerA & markerB) = 0;
            
            VL = circshift(UR,1,2);
            AL = circshift(AR,1,2);
            markerC = ((VL-AL < 0) & (VR-AR > 0));
            markerD = ((VL+AL < 0) & (VR+AR > 0));
            RCOR = markerC.*(VR-AR) + markerD.*(VR+AR);
            LCOR = markerC.*(VL-AL) + markerD.*(VL+AL);
            YCOR = RCOR-LCOR;
            YCOR(markerC & markerD) = 0;
            
            marker = markerA | markerB | markerC | markerD;
        end
        
        function [FX,FY] = AUSMDV_fluxes(FX,FY,U,dx,dy,dt,gamma)
            A = getSoundSpeed(U, gamma);
            W = ToPrimitiveVars(U,gamma);
            
            % x direction
            CM = max(circshift(A,1,1),A);

            PLP = AUSM_SOLVER2D.p_splittingd(W(:,:,end),W(:,:,2),circshift(CM,-1,1),1);
            PRM = AUSM_SOLVER2D.p_splittingd(W(:,:,end),W(:,:,2),CM,-1);
            
            % compute swither s
            K = 10;
            s = min(1,K*abs(PRM-circshift(PLP,1,1))./min(circshift(PLP,1,1),PRM));
            
            % y direction
            CM = max(circshift(A,1,2),A);
            PLP = AUSM_SOLVER2D.p_splittingd(W(:,:,end),W(:,:,3),circshift(CM,-1,2),1);
            PRM = AUSM_SOLVER2D.p_splittingd(W(:,:,end),W(:,:,3),CM,-1);
            s = min(s, min(1,K*abs(PRM-circshift(PLP,1,2))./min(circshift(PLP,1,2),PRM)));
            s = s(2:end,2:end);
            
            % Solution
            [FX,FY] = AUSM_SOLVER2D.AUSMD_fluxes(FX,FY,U,dx,dy,dt,gamma);
            [FX_AUSMV,FY_AUSMV] = AUSM_SOLVER2D.AUSMV_fluxes(FX,FY,U,dx,dy,dt,gamma);
            FX(:,:,2) = 0.5*(1-s).*FX(:,:,2) + 0.5*(1+s).*FX_AUSMV(:,:,2);
            FY(:,:,3) = 0.5*(1-s).*FY(:,:,3) + 0.5*(1+s).*FY_AUSMV(:,:,3);
            
            
            
            % shock-fix
            SHOCKFIX_FLAG = true;
            marker = AUSM_SOLVER2D.IsShockFixLoc(W(:,:,2),W(:,:,3),A);
            S = 0.*A;S(marker) = 1;
            SFcell = S + circshift(S,-1,1) + S + circshift(S,-1,2) ; % x dir sum-up
            SFcell = SFcell * SHOCKFIX_FLAG;
            SFcell = SFcell(2:end,2:end);
            
            [FX_DIS,FY_DIS ]= AUSM_SOLVER2D.AUSM_fluxes(FX,FY,U,dx,dy,dt,gamma);
            FX = (SFcell<=0).*FX + (SFcell>0).*FX_DIS;
            FY = (SFcell<=0).*FY + (SFcell>0).*FY_DIS;
            
            % entropy-fix
            ENTROPYFIX_FLAG = true;
            [marker, XCOR,YCOR] = AUSM_SOLVER2D.IsSonicPoint(W(:,:,2),W(:,:,3),A);
            
            % x direction
            WL = ToPrimitiveVars(circshift(U, 1,1),gamma); %U_{i-1}
            WR = ToPrimitiveVars(U,gamma); % U_i
            
            FL = WL(:,:,1).*compute_psiX(WL,gamma);
            FR = WR(:,:,1).*compute_psiX(WR,gamma);
            FX = FX - 0.125*XCOR(2:end,2:end).*(FR(2:end,2:end,:)-FL(2:end,2:end,:))*ENTROPYFIX_FLAG;
            
            % y direction
            WL = ToPrimitiveVars(circshift(U, 1,2),gamma); %U_{i-1}
            WR = ToPrimitiveVars(U,gamma); % U_i
            
            FL = WL(:,:,1).*compute_psiY(WL,gamma);
            FR = WR(:,:,1).*compute_psiY(WR,gamma);
            FY = FY - 0.125*YCOR(2:end,2:end).*(FR(2:end,2:end,:)-FL(2:end,2:end,:))*ENTROPYFIX_FLAG;
            
        end
        
        function [FX,FY] = MUSCL_AUSM_fluxes(FX,FY,U,dx,dy,dt,gamma)
            A = getSoundSpeed(U, gamma);
            W = ToPrimitiveVars(U,gamma);
            f_limiter = @ospre;
            
            Wm1 = circshift(W, 1,1); %i-1
            Wp1 = circshift(W,-1,1); %i+1
            r1 = (Wp1-W)./(W - Wm1); r1(isnan(r1))=0;
            r2 = (Wm1-W)./(W - Wp1); r2(isnan(r2))=0;
            WLL = W + 0.5*(W - Wm1).*f_limiter(r1);
            WLL = circshift(WLL,1,1);
            WLR = W + 0.5*(W - Wp1).*f_limiter(r2);
            
            AL = getSoundSpeed(ToConservedVars(WLL,gamma), gamma);
            AR = getSoundSpeed(ToConservedVars(WLR,gamma), gamma);
            
            % X direction splitting
            ULP = AUSM_SOLVER2D.u_splitting(WLL(:,:,2),AL, 1);
            URM = AUSM_SOLVER2D.u_splitting(WLR(:,:,2),AR,-1);
            UH = ULP+URM;
            
            PLP = AUSM_SOLVER2D.p_splitting(WLL(:,:,end),WLL(:,:,2),ULP,AL, 1);
            PRM = AUSM_SOLVER2D.p_splitting(WLR(:,:,end),WLR(:,:,2),URM,AR,-1);
            PH = PLP+PRM;
            
            %WL = ToPrimitiveVars(circshift(U, 1,1),gamma); %U_{i-1}
            %WR = ToPrimitiveVars(U,gamma); % U_i
            FL = WLL(:,:,1).*compute_psiX(WLL,gamma);
            FR = WLR(:,:,1).*compute_psiX(WLR,gamma);
            
            FXstar = 0.5*(UH.*(FL+FR) - abs(UH).*(FR-FL));
            FXstar(:,:,2) = FXstar(:,:,2) + PH;
            
            
            % Y direction splitting
            Wm1 = circshift(W, 1,2); %i-1
            Wp1 = circshift(W,-1,2); %i+1
            
            r1 = (Wp1-W)./(W - Wm1); r1(isnan(r1))=0;
            r2 = (Wm1-W)./(W - Wp1); r2(isnan(r2))=0;
            WLL = W + 0.5*(W - Wm1).*f_limiter(r1);
            WLL = circshift(WLL,1,2);
            WLR = W + 0.5*(W - Wp1).*f_limiter(r2);
            
            AL = getSoundSpeed(ToConservedVars(WLL,gamma), gamma);
            AR = getSoundSpeed(ToConservedVars(WLR,gamma), gamma);
            
            ULP = AUSM_SOLVER2D.u_splitting(WLL(:,:,3),AL, 1);
            URM = AUSM_SOLVER2D.u_splitting(WLR(:,:,3),AR,-1);
            UH = ULP+URM;
            
            PLP = AUSM_SOLVER2D.p_splitting(WLL(:,:,end),WLL(:,:,3),ULP,AL, 1);
            PRM = AUSM_SOLVER2D.p_splitting(WLR(:,:,end),WLR(:,:,3),URM,AR,-1);
            PH = PLP+PRM;
            
            %WL = ToPrimitiveVars(circshift(U, 1,2),gamma); %U_{i-1}
            %WR = ToPrimitiveVars(U,gamma); % U_i
            FL = WLL(:,:,1).*compute_psiY(WLL,gamma);
            FR = WLR(:,:,1).*compute_psiY(WLR,gamma);
            
            FYstar = 0.5*(UH.*(FL+FR) - abs(UH).*(FR-FL));
            FYstar(:,:,3) = FYstar(:,:,3) + PH;
            
            
            FX(:,1:end-1,:) = FXstar(2:end,2:end-1,:);
            FY(1:end-1,:,:) = FYstar(2:end-1,2:end,:);
            
            %mustBeReal(flux);
            
        end
        

        function [FX,FY] = AUSM_fluxes(FX,FY,U,dx,dy,dt,gamma)
            A = getSoundSpeed(U, gamma);
            W = ToPrimitiveVars(U,gamma);
            
            % X direction splitting
            ULP = AUSM_SOLVER2D.u_splitting(W(:,:,2),A, 1);
            URM = AUSM_SOLVER2D.u_splitting(W(:,:,2),A,-1);
            UH = circshift(ULP,1,1)+URM;
            
            PLP = AUSM_SOLVER2D.p_splitting(W(:,:,end),W(:,:,2),ULP,A, 1);
            PRM = AUSM_SOLVER2D.p_splitting(W(:,:,end),W(:,:,2),URM,A,-1);
            PH = circshift(PLP,1,1)+PRM;
            
            WL = ToPrimitiveVars(circshift(U, 1,1),gamma); %U_{i-1}
            WR = ToPrimitiveVars(U,gamma); % U_i
            FL = WL(:,:,1).*compute_psiX(WL,gamma);
            FR = WR(:,:,1).*compute_psiX(WR,gamma);

            FXstar = 0.5*(UH.*(FL+FR) - abs(UH).*(FR-FL));
            FXstar(:,:,2) = FXstar(:,:,2) + PH;
            FXstar = FXstar(2:end,:,:);
            
            
            % Y direction splitting
            ULP = AUSM_SOLVER2D.u_splitting(W(:,:,3),A, 1);
            URM = AUSM_SOLVER2D.u_splitting(W(:,:,3),A,-1);
            UH = circshift(ULP,1,2)+URM;
            
            PLP = AUSM_SOLVER2D.p_splitting(W(:,:,end),W(:,:,3),ULP,A, 1);
            PRM = AUSM_SOLVER2D.p_splitting(W(:,:,end),W(:,:,3),URM,A,-1);
            PH = circshift(PLP,1,2)+PRM;
            
            WL = ToPrimitiveVars(circshift(U, 1,2),gamma); %U_{i-1}
            WR = ToPrimitiveVars(U,gamma); % U_i
            FL = WL(:,:,1).*compute_psiY(WL,gamma);
            FR = WR(:,:,1).*compute_psiY(WR,gamma);
            
            FYstar = 0.5*(UH.*(FL+FR) - abs(UH).*(FR-FL));
            FYstar(:,:,3) = FYstar(:,:,3) + PH;
            FYstar = FYstar(:, 2:end,:);
            
            
            FX(:,1:end-1,:) = FXstar(:,2:end-1,:);
            FY(1:end-1,:,:) = FYstar(2:end-1,:,:);
            %FX = FXstar(:,2:end-1,:);
            %FY = FYstar(2:end-1,:,:);
            %mustBeReal(flux);
            
        end
        
        function [FX,FY] = AUSMD_fluxes(FX,FY,U,dx,dy,dt,gamma)
            A = getSoundSpeed(U, gamma);
            W = ToPrimitiveVars(U,gamma);
            
            % X direction
            CM = max(circshift(A,1,1),A);
            POR = W(:,:,end)./W(:,:,1);
            AL = 2*POR./(circshift(POR,-1,1) + POR);
            AR = 2*POR./(circshift(POR, 1,1) + POR);
            
            ULP = AUSM_SOLVER2D.u_splittingd(W(:,:,2),circshift(CM,-1,1),AL, 1);
            URM = AUSM_SOLVER2D.u_splittingd(W(:,:,2),CM,AR,-1);
            UH = circshift(ULP,1,1)+URM;
            
            PLP = AUSM_SOLVER2D.p_splittingd(W(:,:,end),W(:,:,2),circshift(CM,-1,1),1);
            PRM = AUSM_SOLVER2D.p_splittingd(W(:,:,end),W(:,:,2),CM,-1);
            PH = circshift(PLP,1,1)+PRM;
            
            WL = ToPrimitiveVars(circshift(U, 1,1),gamma); %U_{i-1}
            WR = ToPrimitiveVars(U,gamma); % U_i
            RU = circshift(ULP,1,1).*WL(:,:,1) + URM.*WR(:,:,1);
            FL = compute_psiX(WL,gamma);
            FR = compute_psiX(WR,gamma);
            
            FXstar = 0.5*(RU.*(FL+FR) - abs(RU).*(FR-FL));
            FXstar(:,:,2) = FXstar(:,:,2) + PH;
            
            % Y direction
            CM = max(circshift(A,1,2),A);
            POR = W(:,:,end)./W(:,:,1);
            AL = 2*POR./(circshift(POR,-1,2) + POR);
            AR = 2*POR./(circshift(POR, 1,2) + POR);
            
            ULP = AUSM_SOLVER2D.u_splittingd(W(:,:,3),circshift(CM,-1,2),AL, 1);
            URM = AUSM_SOLVER2D.u_splittingd(W(:,:,3),CM,AR,-1);
            UH = circshift(ULP,1,2)+URM;
            
            PLP = AUSM_SOLVER2D.p_splittingd(W(:,:,end),W(:,:,3),circshift(CM,-1,2),1);
            PRM = AUSM_SOLVER2D.p_splittingd(W(:,:,end),W(:,:,3),CM,-1);
            PH = circshift(PLP,1,2)+PRM;
            
            WL = ToPrimitiveVars(circshift(U, 1,2),gamma); %U_{i-1}
            WR = ToPrimitiveVars(U,gamma); % U_i
            RU = circshift(ULP,1,2).*WL(:,:,1) + URM.*WR(:,:,1);
            FL = compute_psiY(WL,gamma);
            FR = compute_psiY(WR,gamma);
            
            FYstar = 0.5*(RU.*(FL+FR) - abs(RU).*(FR-FL));
            FYstar(:,:,3) = FYstar(:,:,3) + PH;
            
            FX(:,1:end-1,:) = FXstar(2:end,2:end-1,:);
            FY(1:end-1,:,:) = FYstar(2:end-1,2:end,:);
            %mustBeReal(flux);
            
        end
        
        function [FX,FY] = AUSMV_fluxes(FX,FY,U,dx,dy,dt,gamma)
            A = getSoundSpeed(U, gamma);
            W = ToPrimitiveVars(U,gamma);
            
            % x direction
            CM = max(circshift(A,1,1),A);
            POR = W(:,:,end)./W(:,:,1);
            AL = 2*POR./(circshift(POR,-1,1) + POR);
            AR = 2*POR./(circshift(POR, 1,1) + POR);
            
            ULP = AUSM_SOLVER2D.u_splittingd(W(:,:,2),circshift(CM,-1,1),AL, 1);
            URM = AUSM_SOLVER2D.u_splittingd(W(:,:,2),CM,AR,-1);
            UH = circshift(ULP,1,1)+URM;
            
            PLP = AUSM_SOLVER2D.p_splittingd(W(:,:,end),W(:,:,2),circshift(CM,-1,1),1);
            PRM = AUSM_SOLVER2D.p_splittingd(W(:,:,end),W(:,:,2),CM,-1);
            PH = circshift(PLP,1,1)+PRM;
            
            WL = ToPrimitiveVars(circshift(U, 1,1),gamma); %U_{i-1}
            WR = ToPrimitiveVars(U,gamma); % U_i
            RU = circshift(ULP,1,1).*WL(:,:,1) + URM.*WR(:,:,1);
            FL = compute_psiX(WL,gamma);
            FR = compute_psiX(WR,gamma);
            
            FXstar = 0.5*(RU.*(FL+FR) - abs(RU).*(FR-FL));
            FXstar(:,:,2) = circshift(ULP,1,1).*WL(:,:,1).*WL(:,:,2) + URM.*WR(:,:,1).*WR(:,:,2) + PH;
            
            % y direction
            CM = max(circshift(A,1,2),A);
            POR = W(:,:,end)./W(:,:,1);
            AL = 2*POR./(circshift(POR,-1,2) + POR);
            AR = 2*POR./(circshift(POR, 1,2) + POR);
            
            ULP = AUSM_SOLVER2D.u_splittingd(W(:,:,3),circshift(CM,-1,2),AL, 1);
            URM = AUSM_SOLVER2D.u_splittingd(W(:,:,3),CM,AR,-1);
            UH = circshift(ULP,1,2)+URM;
            
            PLP = AUSM_SOLVER2D.p_splittingd(W(:,:,end),W(:,:,3),circshift(CM,-1,2),1);
            PRM = AUSM_SOLVER2D.p_splittingd(W(:,:,end),W(:,:,3),CM,-1);
            PH = circshift(PLP,1,2)+PRM;
            
            WL = ToPrimitiveVars(circshift(U, 1,2),gamma); %U_{i-1}
            WR = ToPrimitiveVars(U,gamma); % U_i
            RU = circshift(ULP,1,2).*WL(:,:,1) + URM.*WR(:,:,1);
            FL = compute_psiY(WL,gamma);
            FR = compute_psiY(WR,gamma);
            
            FYstar = 0.5*(RU.*(FL+FR) - abs(RU).*(FR-FL));
            FYstar(:,:,3) = circshift(ULP,1,2).*WL(:,:,1).*WL(:,:,3) + URM.*WR(:,:,1).*WR(:,:,3) + PH;
            
            FX(:,1:end-1,:) = FXstar(2:end,2:end-1,:);
            FY(1:end-1,:,:) = FYstar(2:end-1,2:end,:);
            %mustBeReal(flux);
            
        end
        
        
        
    end
end

%% limiter
function phi = minmod(r)
phi = max(0,min(1,r));
end
function phi = superbee(r)
phi = max(0,max(min(1,2*r),min(r,2)));
end
function phi = mc(r)
phi = max(0,min(2,min(2*r,0.5*(1+r))));
end

function phi = ospre(r)
phi = 1.5*(r.^2+r)./(r.^2+r+1);
phi(isinf(r)) = 1.5;
end

function phi = vanleer(r)
phi = (r+abs(r))/(1+abs(r));
phi(isinf(r)) = 2;
end
%% Helper functions
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
W = 0.*U;
W(:,:,1) = U(:,:,1);
W(:,:,2:(end-1)) = U(:,:,2:(end-1))./U(:,:,1);
W(:,:,end) = (U(:,:,end) - 0.5.*(getVelMagSqr(U)).*U(:,:,1))*(gamma-1);
end

function U = ToConservedVars(W, gamma)
U = 0.*W;
U(:,:,1) = W(:,:,1);
U(:,:,2:(end-1)) = W(:,:,1).*W(:,:,2:(end-1));
U(:,:,end) = 0.5*W(:,:,1).*(sum(W(:,:,2:(end-1)).^2,3)) + W(:,:,end)/(gamma-1);
end

function FU = compute_psiX(W,gamma)
FU = cat(3,  1+0.*W(:,:,1) ...
    ,W(:,:,2) ...
    ,W(:,:,3)...
    ,(W(:,:,4)./W(:,:,1) + (0.5*(W(:,:,2).^2+W(:,:,3).^2) + W(:,:,4)./W(:,:,1)/(gamma-1))) );

end

function FV = compute_psiY(W,gamma)
FV = cat(3,  1+0.*W(:,:,1) ...
    ,W(:,:,2)...
    ,W(:,:,3)...
    ,(W(:,:,4)./W(:,:,1) + (0.5*(W(:,:,2).^2+W(:,:,3).^2) + W(:,:,4)./W(:,:,1)/(gamma-1))) );
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