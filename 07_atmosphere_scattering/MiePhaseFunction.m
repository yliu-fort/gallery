function out = MiePhaseFunction(g, nu) 
  k = 3.0 ./ (8.0 * pi) .* (1.0 - g .* g) ./ (2.0 + g .* g);
  out = k .* (1.0 + nu .* nu) ./ (1.0 + g .* g - 2.0 * g .* nu) .^ 1.5;
end