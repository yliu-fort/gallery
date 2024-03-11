function out = DistanceToTopAtmosphereBoundary(atmosphere, r, mu)
  %assert(all(r(:) <= atmosphere.Rt));
  %assert(all(mu(:) >= -1.0 & mu(:) <= 1.0));
discriminant = r.^2.*(mu.^2-1) + atmosphere.Rt^2;
out = ClampDistance(-r.*mu + sqrt(discriminant));
end