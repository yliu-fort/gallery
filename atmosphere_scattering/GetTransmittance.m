function out = GetTransmittance(...
    atmosphere,...
    transmittance_texture,...
    r, mu, d, ray_r_mu_intersects_ground) 
  %assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
  %assert(mu >= -1.0 && mu <= 1.0);
  %assert(d >= 0.0 * m);
  mask = ray_r_mu_intersects_ground;
  r_d = ClampRadius(sqrt(d .* d + 2.0 * r .* mu .* d + r .* r), atmosphere.Rg, atmosphere.Rt);
  mu_d = ClampCosine((r .* mu + d) ./ r_d);
  [u_r, u_mu] = RMuEncoding(atmosphere, r, mu);
  [u_rd, u_mud] = RMuEncoding(atmosphere, r_d, mu_d);
  [~, u_mu1] = RMuEncoding(atmosphere, r, -mu);
  [~, u_mud1] = RMuEncoding(atmosphere, r_d, -mu_d);

    out(mask,:) = min(...
        transmittance_texture(u_rd(mask), u_mud1(mask)) ./...
        transmittance_texture(u_r(mask), u_mu1(mask)),...
        1.0);
    out(~mask,:) = min(...
        transmittance_texture(u_r(~mask), u_mu(~mask)) ./...
        transmittance_texture(u_rd(~mask), u_mud(~mask)),...
        1.0);
end