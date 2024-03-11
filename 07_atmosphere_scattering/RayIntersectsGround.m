function out = RayIntersectsGround(r, mu, Rg)
  out = mu < 0.0 & r .* r .* (mu .* mu - 1.0) + ...
      Rg .* Rg >= 0.0;
end