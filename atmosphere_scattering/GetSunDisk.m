function out = GetSunDisk(atmosphere,r, mu_s) 
  sin_theta_h = atmosphere.Rg ./ r;
  cos_theta_h = -sqrt(max(1.0 - sin_theta_h .* sin_theta_h, 0.0));

    out = ...
      smoothstep(-sin_theta_h * atmosphere.sun_angular_radius, ...
                 sin_theta_h * atmosphere.sun_angular_radius, ...
                 mu_s - cos_theta_h);
end