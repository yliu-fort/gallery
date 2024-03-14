classdef AUSM_SOLVER1D
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
        
        function marker = IsShockFixLoc(UR,AR)
            UL = circshift(UR,1,1);
            AL = circshift(AR,1,1);
            marker = ((UL-AL > 0) & (UR-AR < 0)) | ((UL+AL > 0) & (UR+AR < 0));
            marker = marker & circshift(marker,1,1);
        end
        
        function [marker, RCOR,LCOR] = IsSonicPoint(UR,AR)
            UL = circshift(UR,1,1);
            AL = circshift(AR,1,1);
            markerA = ((UL-AL < 0) & (UR-AR > 0));
            markerB = ((UL+AL < 0) & (UR+AR > 0));
            marker = 1.*markerA - 1.*markerB;
            RCOR = markerA.*(UR-AR) + markerB.*(UR+AR);
            LCOR = markerA.*(UL-AL) + markerB.*(UL+AL);
            RCOR(markerA & markerB) = 0;
            LCOR(markerA & markerB) = 0;
        end
        
        function flux = AUSMDV_fluxes(flux,U,dx,dt,gamma,nvars)
            A = getSoundSpeed(U, gamma);
            W = ToPrimitiveVars(U,gamma);
            
            CM = max(circshift(A,1,1),A);
            
            PLP = AUSM_SOLVER1D.p_splittingd(W(:,3),W(:,2),circshift(CM,-1,1),1);
            PRM = AUSM_SOLVER1D.p_splittingd(W(:,3),W(:,2),CM,-1);
            
            % compute swither s
            K = 10;
            s = min(1,K*abs(PRM-circshift(PLP,1,1))./min(circshift(PLP,1,1),PRM));
            %s(~isfinite(s))=0;
            s = s(2:end,:);
            
            WL = ToPrimitiveVars(circshift(U, 1,1),gamma); %U_{i-1}
            WR = ToPrimitiveVars(U,gamma); % U_i
            
            FL = [WL(:,1) WL(:,1).*WL(:,2) (WL(:,3) + WL(:,1).*0.5.*WL(:,2).^2 + WL(:,3)/(gamma-1))];
            FR = [WR(:,1) WR(:,1).*WR(:,2) (WR(:,3) + WR(:,1).*0.5.*WR(:,2).^2 + WR(:,3)/(gamma-1))];
            
            % Solution
            flux = AUSM_SOLVER1D.AUSMD_fluxes(flux,U,dx,dt,gamma,nvars);
            flux_AUSMV = AUSM_SOLVER1D.AUSMV_fluxes(flux,U,dx,dt,gamma,nvars);
            flux(:,2) = 0.5*(1-s).*flux(:,2) + 0.5*(1+s).*flux_AUSMV(:,2);
            
            % shock-fix
            SHOCKFIX_FLAG = true;
            marker = AUSM_SOLVER1D.IsShockFixLoc(W(:,2),A);
            S = 0.*A;S(marker) = 1;
            SFcell = S + circshift(S,-1,1); % x dir sum-up
            SFcell = SFcell * SHOCKFIX_FLAG;
            SFcell = SFcell(2:end,:);
            
            DISflux = AUSM_SOLVER1D.AUSM_fluxes(flux,U,dx,dt,gamma,nvars);
            flux(SFcell>0,:) = DISflux(SFcell>0,:);
            
            % entropy-fix
            ENTROPYFIX_FLAG = true;
            [marker, RCOR,LCOR] = AUSM_SOLVER1D.IsSonicPoint(W(:,2),A);
            flux = flux - 0.125*(RCOR(2:end,:)-LCOR(2:end,:)).*(FR(2:end,:)-FL(2:end,:))*ENTROPYFIX_FLAG;
            
        end
        
        function flux = AUSM_fluxes(flux,U,dx,dt,gamma,nvars)
            A = getSoundSpeed(U, gamma);
            W = ToPrimitiveVars(U,gamma);
            
            ULP = AUSM_SOLVER1D.u_splitting(W(:,2),A, 1);
            URM = AUSM_SOLVER1D.u_splitting(W(:,2),A,-1);
            UH = circshift(ULP,1,1)+URM;
            
            PLP = AUSM_SOLVER1D.p_splitting(W(:,3),W(:,2),ULP,A, 1);
            PRM = AUSM_SOLVER1D.p_splitting(W(:,3),W(:,2),URM,A,-1);
            PH = circshift(PLP,1,1)+PRM;
            
            WL = ToPrimitiveVars(circshift(U, 1,1),gamma); %U_{i-1}
            WR = ToPrimitiveVars(U,gamma); % U_i
            %FL = [WL(:,1) WL(:,1).*WL(:,2) (WL(:,3) + WL(:,1).*0.5.*WL(:,2).^2 + WL(:,3)/(gamma-1))];
            %FR = [WR(:,1) WR(:,1).*WR(:,2) (WR(:,3) + WR(:,1).*0.5.*WR(:,2).^2 + WR(:,3)/(gamma-1))];
            FL = [WL(:,1) WL(:,1).*WL(:,2) (WL(:,3) + WL(:,1).*0.5.*WL(:,2).^2 + WL(:,3)/(gamma-1))];
            FR = [WR(:,1) WR(:,1).*WR(:,2) (WR(:,3) + WR(:,1).*0.5.*WR(:,2).^2 + WR(:,3)/(gamma-1))];
            
            flux = 0.5*(UH.*(FL+FR) - abs(UH).*(FR-FL));
            flux(:,2) = flux(:,2) + PH;
            flux = flux(2:end,:);
            
            %mustBeReal(flux);
            
        end
        
        function flux = AUSMD_fluxes(flux,U,dx,dt,gamma,nvars)
            A = getSoundSpeed(U, gamma);
            W = ToPrimitiveVars(U,gamma);
            
            CM = max(circshift(A,1,1),A);
            POR = W(:,3)./W(:,1);
            AL = 2*POR./(circshift(POR,-1,1) + POR);
            AR = 2*POR./(circshift(POR, 1,1) + POR);
            
            ULP = AUSM_SOLVER1D.u_splittingd(W(:,2),circshift(CM,-1,1),AL, 1);
            URM = AUSM_SOLVER1D.u_splittingd(W(:,2),CM,AR,-1);
            UH = circshift(ULP,1,1)+URM;
            
            PLP = AUSM_SOLVER1D.p_splittingd(W(:,3),W(:,2),circshift(CM,-1,1),1);
            PRM = AUSM_SOLVER1D.p_splittingd(W(:,3),W(:,2),CM,-1);
            PH = circshift(PLP,1,1)+PRM;
            
            WL = ToPrimitiveVars(circshift(U, 1,1),gamma); %U_{i-1}
            WR = ToPrimitiveVars(U,gamma); % U_i
            RU = circshift(ULP,1,1).*WL(:,1) + URM.*WR(:,1);
            FL = [1+0*WL(:,1) WL(:,2) 0.5*WL(:,2).^2 + WL(:,3)./WL(:,1)*gamma/(gamma-1)];
            FR = [1+0*WR(:,1) WR(:,2) 0.5*WR(:,2).^2 + WR(:,3)./WR(:,1)*gamma/(gamma-1)];
            
            flux = 0.5*(RU.*(FL+FR) - abs(RU).*(FR-FL));
            flux(:,2) = flux(:,2) + PH;
            flux = flux(2:end,:);
            
            mustBeReal(flux);
            
        end
        
        function flux = AUSMV_fluxes(flux,U,dx,dt,gamma,nvars)
            A = getSoundSpeed(U, gamma);
            W = ToPrimitiveVars(U,gamma);
            
            CM = max(circshift(A,1,1),A);
            POR = W(:,3)./W(:,1);
            AL = 2*POR./(circshift(POR,-1,1) + POR);
            AR = 2*POR./(circshift(POR, 1,1) + POR);
            
            ULP = AUSM_SOLVER1D.u_splittingd(W(:,2),circshift(CM,-1,1),AL, 1);
            URM = AUSM_SOLVER1D.u_splittingd(W(:,2),CM,AR,-1);
            UH = circshift(ULP,1,1)+URM;
            
            PLP = AUSM_SOLVER1D.p_splittingd(W(:,3),W(:,2),circshift(CM,-1,1),1);
            PRM = AUSM_SOLVER1D.p_splittingd(W(:,3),W(:,2),CM,-1);
            PH = circshift(PLP,1,1)+PRM;
            
            WL = ToPrimitiveVars(circshift(U, 1,1),gamma); %U_{i-1}
            WR = ToPrimitiveVars(U,gamma); % U_i
            RU = circshift(ULP,1,1).*WL(:,1) + URM.*WR(:,1);
            FL = [1+0*WL(:,1) WL(:,2) 0.5*WL(:,2).^2 + WL(:,3)./WL(:,1)*gamma/(gamma-1)];
            FR = [1+0*WR(:,1) WR(:,2) 0.5*WR(:,2).^2 + WR(:,3)./WR(:,1)*gamma/(gamma-1)];
            
            flux = 0.5*(RU.*(FL+FR) - abs(RU).*(FR-FL));
            flux(:,2) = circshift(ULP,1,1).*WL(:,1).*WL(:,2) + URM.*WR(:,1).*WR(:,2) + PH;
            flux = flux(2:end,:);
            
            mustBeReal(flux);
            
        end
        

        
    end
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
