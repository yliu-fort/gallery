function [u_r, u_mu, u_mu_s, u_nu] = RMuMusNuEncoding(atmosphere, r, mu, mu_s, nu, ray_r_mu_intersects_ground)
    % Distance to top atmosphere boundary for a horizontal ray at ground level.
    Rg = atmosphere.Rg;
    Rt = atmosphere.Rt;
    H = sqrt(Rt^2 - Rg^2);
    %mu_s_min = -sqrt(1-(Rg/Rt)^2);
    mu_s_min = atmosphere.mu_s_min;
    % Distance to the horizon.
    rho = sqrt(r.^2 - Rg^2);
    u_r = rho ./ H;

    % Discriminant of the quadratic equation for the intersections of the ray
    % (r,mu) with the ground (see RayIntersectsGround).
    r_mu = r .* mu;
    discriminant = r_mu.^2 - r.^2 + Rg^2;
    if nargin < 7
        ray_r_mu_intersects_ground = RayIntersectsGround(r, mu, Rg);
    end

    % Distance to the ground for the ray (r,mu), and its minimum and maximum
    % values over all mu - obtained for (r,-1) and (r,mu_horizon).
    d = -r_mu - SafeSqrt(discriminant);
    d_min = r - Rg;
    d_max = rho;
    u_mu1 = 0.5 - 0.5 * (d - d_min) ./ (d_max - d_min);
    u_mu1(d_max == d_min) = 0;

    % Distance to the top atmosphere boundary for the ray (r,mu), and its
    % minimum and maximum values over all mu - obtained for (r,1) and
    % (r,mu_horizon).
    d = -r_mu + SafeSqrt(discriminant + H^2);
    d_min = Rt - r;
    d_max = rho + H;
    u_mu = 0.5 + 0.5 * (d - d_min) ./ (d_max - d_min);
    u_mu(ray_r_mu_intersects_ground) = u_mu1(ray_r_mu_intersects_ground);

    d = DistanceToTopAtmosphereBoundary(atmosphere, Rg, mu_s);
    d_min = Rt - Rg;
    d_max = H;
    a = (d - d_min) ./ (d_max - d_min);
    D = DistanceToTopAtmosphereBoundary(atmosphere, Rg, mu_s_min);
    A = (D - d_min) ./ (d_max - d_min);
    % An ad-hoc function equal to 0 for mu_s = mu_s_min (because then d = D and
    % thus a = A), equal to 1 for mu_s = 1 (because then d = d_min and thus
    % a = 0), and with a large slope around mu_s = 0, to get more texture 
    % samples near the horizon.
    u_mu_s = max(1.0 - a / A, 0.0) ./ (1.0 + a);

    u_nu = (nu + 1.0) / 2.0;

end