function [rayleigh, mie] = ComputeSingleScattering(...
    atmosphere,...
    r, mu, mu_s, nu)

ray_r_mu_intersects_ground = RayIntersectsGround(r, mu, atmosphere.Rg);

% Number of intervals for the numerical integration.
SAMPLE_COUNT = 50;
% The integration step, i.e. the length of each integration interval.
dx = arrayfun(@(r, mu, ray_r_mu_intersects_ground)...
    DistanceToNearestAtmosphereBoundary(atmosphere, r, mu, ray_r_mu_intersects_ground),...
    r, mu,ray_r_mu_intersects_ground) ./ SAMPLE_COUNT;

rayleigh_sum = 0;
mie_sum = 0;
for i = 0:SAMPLE_COUNT
    d_i = i * dx;
    r_d = ClampRadius(sqrt(d_i .^ 2 + 2.0 * r .* mu .* d_i + r .^ 2), atmosphere.Rg, atmosphere.Rt);
    mu_s_d = ClampCosine((r .* mu_s + d_i .* nu) ./ r_d);

    transmittence = ...
      GetTransmittance( ...
          atmosphere, atmosphere.cachedTransmittence, r, mu, d_i, ray_r_mu_intersects_ground) .* ...
      GetTransmittanceToSun( ...
          atmosphere, atmosphere.cachedTransmittence, r_d, mu_s_d);

    weight_i = 1.0;
    if (i == 0 || i == SAMPLE_COUNT), weight_i = 0.5; end
    rayleigh_sum = rayleigh_sum + weight_i * atmosphere.rayleigh_density(r_d - atmosphere.Rg) .* transmittence;
    mie_sum = mie_sum + weight_i * atmosphere.mie_density(r_d - atmosphere.Rg) .* transmittence;
end
rayleigh = rayleigh_sum .* atmosphere.rayleigh_scattering.' .* dx;
mie = mie_sum .* atmosphere.mie_scattering.' .* dx;
end