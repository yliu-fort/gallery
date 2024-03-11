function nu = ClampNu(mu, mu_s, nu)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nu = max(min(nu, mu .* mu_s - sqrt((1.0 - mu .* mu) .* (1.0 - mu_s .* mu_s))),...
      mu .* mu_s + sqrt((1.0 - mu .* mu) .* (1.0 - mu_s .* mu_s)));
end