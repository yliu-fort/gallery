clear;clc;close all;

cache_name = "Earth";
output_path = "../pysr_core/";

% Stage 0 Solar physics
kConstantSolarIrradiance = 1.5;
%kLambda = [680;550;440]*1e-3; % wavelength, *1e3 nm
kLambda = [360:10:830].'*1e-3;
kSolarIrradiance = 0*kLambda + kConstantSolarIrradiance;
kSunAngularRadius = 0.00935/2;
kSolarRadiance = kSolarIrradiance ./ (pi * kSunAngularRadius.^2);

% Stage 1: Scattering physics
switch cache_name
    case 'Mars'
        CustomScaleCoeff = 10.0;
        Rscale = 3389.5e3; % 3397.5 km
        AtmosphereHeight = 110e3*CustomScaleCoeff; % 119 km
        kMiePhaseFunctionG = [0.1;0.2;0.8];
        N = 2.5e25*0.02 / sqrt(CustomScaleCoeff); % molecule density at sea level/ mol
        H0 = 10800/Rscale*CustomScaleCoeff; % Atmosphere thickness, m
        H1 = 5500/Rscale*CustomScaleCoeff;
        kAtmosphereGasFractions = [0.0,0.03,0.02,0,0,0.95,0];
    case 'Earth'
        CustomScaleCoeff = 10.0;
        Rscale = 6371e3; % 6371 km
        AtmosphereHeight = 85e3*CustomScaleCoeff; % 85 km
        kMiePhaseFunctionG = [0.8;0.8;0.8];
        N = 2.5e25 / sqrt(CustomScaleCoeff); % molecule density at sea level/ mol
        H0 = 8500/Rscale*CustomScaleCoeff; % Atmosphere thickness, m
        H1 = 1100/Rscale*CustomScaleCoeff;
        kAtmosphereGasFractions = [0.209,0.781,0.0093,0,0,0.0004,0];
end

Ns = 6.02e23; % Avogadro's number
kRayleigh = (N/Ns) * getCrossSectionTypicalGases(kAtmosphereGasFractions, kLambda*1e3) * Rscale;

switch cache_name
    case 'Mars'
    kMieAngstromAlpha = -2.85; % -2.85 for mars
    case 'Earth'
    kMieAngstromAlpha = 0.0;
end
kMieAngstromBeta = 5.328e-3 / 41.5282 * (N/Ns);
kMieSingleScatteringAlbedo = 0.9;

kMie = kMieAngstromBeta ./ H1 .* ((kLambda*1e3) .^ -kMieAngstromAlpha);
switch cache_name
    case 'Mars'
    kMieUpperBoundCoeff = 20.0; % ratio of dust scattering intensity / air molecule intensity
    kMieUpperBound = kMieUpperBoundCoeff * 28.2872 * (N/Ns) / 41.5282;
    kMie = kMie / max(kMie) * kMieUpperBound;
end
% Stage 2: Space directions
planetCenter = [0,0,0];
rGround = Rscale; % scaled to 1
rAtmosphere = Rscale + AtmosphereHeight; % scaled to 0.013 for earth

% Stage 3: Setup atmosphere datastructure
atmosphere.Rg = rGround / Rscale;
atmosphere.Rt = rAtmosphere / Rscale;
atmosphere.sun_angular_radius = kSunAngularRadius;
atmosphere.rayleigh_density = @(h)exp(-(h)/H0);
atmosphere.rayleigh = kRayleigh;
atmosphere.rayleigh_scattering = kRayleigh ./ kLambda.^4;
atmosphere.mie_density = @(h)exp(-(h)/H1);
atmosphere.mie_scattering = kMie * kMieSingleScatteringAlbedo;
atmosphere.mie_extinction = kMie;
atmosphere.mie_phase_function_g = kMiePhaseFunctionG;
atmosphere.mu_s_min = -sqrt(1-(rGround/rAtmosphere)^2);

%% Stage 3.1: Precompute optical depth
u_rDist = linspace(0,1,64);
u_muDist = linspace(0,1,256);
[U_R0,U_MU0] = ndgrid(u_rDist, u_muDist);
[R,MU] = RMuDecoding(atmosphere, U_R0,U_MU0);

for i = 1:size(R,1)
    for j = 1:size(R,2)
        Transmittence = ComputeTransmittanceToTopAtmosphereBoundary(atmosphere, R(i,j), MU(i,j));
        for l = 1:numel(kLambda)
           Tex_Transmittence(i,j,l) = Transmittence(l);
        end
    end
end

atmosphere.cachedTransmittence = @(u_r, u_mu)(interpn(U_R0,U_MU0,Tex_Transmittence,u_r,u_mu,'linear',0));

%% 3.2 Precompute transmittence

u_rDist = linspace(0,1,32);
u_muDist = [linspace(1e-6,0.5-1e-6,64),linspace(0.5+1e-6,1-1e-6,64)];
u_musDist = linspace(0,1,32);
u_nuDist = linspace(0,1,8);
[U_R,U_MU,U_MUS,U_NU] = ndgrid(u_rDist,u_muDist,u_musDist,u_nuDist);
[R,MU,MUS,NU, ~] = RMuMusNuDecoding(atmosphere, U_R, U_MU, U_MUS, U_NU);

% transmittence
[rayleigh, mie] = ComputeSingleScattering(atmosphere, ...
    R(:), MU(:), MUS(:), NU(:));
for l = 1:numel(kLambda)
    Tex_Rayleigh(:,:,:,:,l) = reshape(rayleigh(:,l), size(R));
    Tex_Mie(:,:,:,:,l) = reshape(mie(:,l), size(R));
end

%% Integrate to table
TR = [];
SCATTER = [];
for l = 1:numel(kLambda)
    TR = [TR;[reshape(U_R0,[],1) ...
        reshape(U_MU0,[],1) ...
        kLambda(l) + zeros(numel(U_R0),1) ...
        reshape(Tex_Transmittence(:,:,l),[],1)]];
    SCATTER = [SCATTER;[reshape(U_R,[],1) ...
        reshape(U_MU,[],1) ...
        reshape(U_MUS,[],1) ...
        reshape(U_NU,[],1) ...
        kLambda(l) + zeros(numel(U_R),1) ...
        reshape(Tex_Rayleigh(:,:,:,:,l),[],1) ...
        reshape(Tex_Mie(:,:,:,:,l),[],1)]];
end

save(output_path+cache_name+"_tr.mat", "TR")
save(output_path+cache_name+".mat", "SCATTER", "-V7.3")

%% print hdr
Tex_TransmittenceC = Tex_Transmittence(:,:,[33,20,9]);
Tex_RayleighC = permute(Tex_Rayleigh(:,:,:,:,[33,20,9]), [2,1,3,4,5]);
Tex_MieC = permute(Tex_Mie(:,:,:,:,[33,20,9]), [2,1,3,4,5]);
Tex_Rayleigh_C(:,:,1) = reshape(Tex_RayleighC(:,:,:,:,1),32*128,[]);
Tex_Rayleigh_C(:,:,2) = reshape(Tex_RayleighC(:,:,:,:,2),32*128,[]);
Tex_Rayleigh_C(:,:,3) = reshape(Tex_RayleighC(:,:,:,:,3),32*128,[]);
Tex_Mie_C(:,:,1) = reshape(Tex_MieC(:,:,:,:,1),32*128,[]);
Tex_Mie_C(:,:,2) = reshape(Tex_MieC(:,:,:,:,2),32*128,[]);
Tex_Mie_C(:,:,3) = reshape(Tex_MieC(:,:,:,:,3),32*128,[]);
hdrwrite(Tex_TransmittenceC,output_path+cache_name+'_tr.hdr')
hdrwrite(Tex_Mie_C,output_path+cache_name+'_mie.hdr')
hdrwrite(Tex_Rayleigh_C,output_path+cache_name+'_rayleigh.hdr')
