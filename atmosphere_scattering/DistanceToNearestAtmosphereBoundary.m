function out = DistanceToNearestAtmosphereBoundary(atmosphere, r, mu, ray_r_mu_intersects_ground) 
  if (ray_r_mu_intersects_ground) 
    out = DistanceToBottomAtmosphereBoundary(atmosphere, r, mu);
   else 
    out = DistanceToTopAtmosphereBoundary(atmosphere, r, mu);
  end
end