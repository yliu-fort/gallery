function [r, mu, mu_s, nu, ray_r_mu_intersects_ground] = RMuMusNuDecoding(atmosphere, u_r, u_mu, u_mu_s, u_nu)
    % Distance to top atmosphere boundary for a horizontal ray at ground level.
    Rg = atmosphere.Rg;
    Rt = atmosphere.Rt; 
    H = sqrt(Rt^2 - Rg^2);
    %mu_s_min = -sqrt(1-(Rg/Rt)^2);
    mu_s_min = atmosphere.mu_s_min;
    % Distance to the horizon.
    rho = H .* u_r;
    r = sqrt(rho.^2 + Rg^2);
%{
      if (uvwz.z < 0.5) {
    // Distance to the ground for the ray (r,mu), and its minimum and maximum
    // values over all mu - obtained for (r,-1) and (r,mu_horizon) - from which
    // we can recover mu:
    Length d_min = r - atmosphere.bottom_radius;
    Length d_max = rho;
    Length d = d_min + (d_max - d_min) * GetUnitRangeFromTextureCoord(
        1.0 - 2.0 * uvwz.z, SCATTERING_TEXTURE_MU_SIZE / 2);
    mu = d == 0.0 * m ? Number(-1.0) :
        ClampCosine(-(rho * rho + d * d) / (2.0 * r * d));
    ray_r_mu_intersects_ground = true;
  } else {
    // Distance to the top atmosphere boundary for the ray (r,mu), and its
    // minimum and maximum values over all mu - obtained for (r,1) and
    // (r,mu_horizon) - from which we can recover mu:
    Length d_min = atmosphere.top_radius - r;
    Length d_max = rho + H;
    Length d = d_min + (d_max - d_min) * GetUnitRangeFromTextureCoord(
        2.0 * uvwz.z - 1.0, SCATTERING_TEXTURE_MU_SIZE / 2);
    mu = d == 0.0 * m ? Number(1.0) :
        ClampCosine((H * H - rho * rho - d * d) / (2.0 * r * d));
    ray_r_mu_intersects_ground = false;
  }
%}
    % Distance to the ground for the ray (r,mu), and its minimum and maximum
    % values over all mu - obtained for (r,-1) and (r,mu_horizon) - from which
    % we can recover mu:
    % (Todo)what happens when r = 0 needs to be inspected.
    d_min = r - Rg;
    d_max = rho;
    d = d_min + (d_max - d_min) .* (1.0 - 2.0 * u_mu);
    mu1 = ClampCosine(-(rho.^2 + d.^2) ./ (2.0 * r .* d));
    mu1(d == 0) = -1;

    % Distance to the top atmosphere boundary for the ray (r,mu), and its
    % minimum and maximum values over all mu - obtained for (r,1) and
    % (r,mu_horizon) - from which we can recover mu:
    d_min = Rt - r;
    d_max = rho + H;
    d = d_min + (d_max - d_min) .* (2.0 * u_mu - 1.0);
    mu = ClampCosine((H.^2 - rho.^2 - d.^2) ./ (2.0 * r .* d));
    mu(d == 0) = 1;
    
    mu(u_mu < 0.5) = mu1(u_mu < 0.5);
    ray_r_mu_intersects_ground = u_mu < 0.5;


    d_min = Rt - Rg;
    d_max = H;
    D = DistanceToTopAtmosphereBoundary(atmosphere, Rg, mu_s_min);
    A = (D - d_min) ./ (d_max - d_min);
    a = (A - u_mu_s .* A) ./ (1.0 + u_mu_s .* A);
    d = d_min + min(a, A) .* (d_max - d_min);
    mu_s = ClampCosine((H.^2 - d.^2) ./ (2.0 * Rg .* d));
    mu_s(d == 0) = 1;
    nu = ClampCosine(u_nu * 2.0 - 1.0);
end
