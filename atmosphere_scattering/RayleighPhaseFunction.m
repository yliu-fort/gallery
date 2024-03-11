function out = RayleighPhaseFunction(nu) 
  k = 3.0 ./ (16.0 * pi);
  out = k .* (1.0 + nu .* nu);
end