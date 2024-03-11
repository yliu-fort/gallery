function out = DistanceToBottomAtmosphereBoundary(atmosphere, ...
    r, mu) 
  assert(all(r >= atmosphere.Rg));
  assert(all(mu >= -1.0 & mu <= 1.0));
  discriminant = r .* r .* (mu .* mu - 1.0) + ...
      atmosphere.Rg .* atmosphere.Rg;
  out = ClampDistance(-r .* mu - SafeSqrt(discriminant));
end