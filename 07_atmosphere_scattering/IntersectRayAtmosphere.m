function [intersections, hits_ground] = IntersectRayAtmosphere(viewOrigin, ...
                                    viewDirection, ...
                                    sphereCenters, ...
                                    sphereRadius) 
% return a series of intersections to the spheres if exists.
% the ray ends if collides with the ground or out of the atmosphere.
nSpheres = size(sphereCenters,1);
intersections = [];
hits_ground = false;
tolOnSphere = 1e-10;

% find nearest collision to the atmos layer
% if collide with ground or nothing, end the loop.
while 1
    distance = nan;
    intersection = zeros(1,3)+nan;
    hitIndex = 0;
    for i = 1:nSpheres
        [t, q] = IntersectRaySphere(viewOrigin, ...
                                        viewDirection, ...
                                        sphereCenters(i,:), ...
                                        sphereRadius(i));
        if (t > tolOnSphere) && (isnan(distance) || (t < distance))
            distance = t;
            intersection = q;
            hitIndex = i;
        end
    end

    % return empty if hit nothing
    if isnan(distance), break; end

    intersections = [intersections;intersection];
    viewOrigin = intersection;

    % return if hit ground
    if hitIndex == 1
        hits_ground = true;
        break;
    end

end

end