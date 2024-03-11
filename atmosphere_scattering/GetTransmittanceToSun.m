function out = GetTransmittanceToSun(atmosphere,transmittance_texture,r, mu_s) 
  %[u_r, u_mus] = RMuEncoding(atmosphere, r, mu_s);
  %out(~hit_ground, :) = transmittance_texture(u_r(~hit_ground), u_mus(~hit_ground));

  out = squeeze(GetTransmittanceToTopAtmosphereBoundary( ...
          atmosphere, transmittance_texture, r, mu_s)) .* ...
      GetSunDisk(atmosphere, r, mu_s);
end