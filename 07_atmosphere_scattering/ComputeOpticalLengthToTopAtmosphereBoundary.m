function result = ComputeOpticalLengthToTopAtmosphereBoundary(...
    atmosphere, profile, r, mu) 
  % Number of intervals for the numerical integration.
  SAMPLE_COUNT = 500;
  % The integration step, i.e. the length of each integration interval.
  dx = DistanceToTopAtmosphereBoundary(atmosphere, r, mu) / SAMPLE_COUNT;
  % Integration loop.
  result = 0.0; % m
  for i = 0:SAMPLE_COUNT
    d_i = i .* dx;
    % Distance between the current sample point and the planet center.
    r_i = sqrt(d_i .* d_i + 2.0 .* r .* mu .* d_i + r .* r);
    % Number density at the current sample point (divided by the number density
    % at the bottom of the atmosphere, yielding a dimensionless number).
    y_i = profile(r_i - atmosphere.Rg);
    % Sample weight (from the trapezoidal rule).
    weight_i = 1.0;
    if (i == 0 || i == SAMPLE_COUNT), weight_i = 0.5; end
    result = result + y_i .* weight_i .* dx;
  end
end