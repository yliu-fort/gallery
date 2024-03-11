function [distance, intersection] = IntersectRaySphere(viewOrigin, ...
                                    viewDirection, ...
                                    sphereCenter, ...
                                    sphereRadius) 
nPixels = 1;

b_inside_circle = sqrt(sum((viewOrigin - sphereCenter).^2,2)) <= sphereRadius;

a = dot(viewDirection, viewDirection, 2);

b = 2 * dot(viewDirection, (viewOrigin - sphereCenter),2);

c = sum((viewOrigin - sphereCenter).^2,2) - sphereRadius^2;

discriminant = b.^2 - 4 .* a .* c;

% if the ray hits the sphere
hits = discriminant >= 0;

distance1 = NaN(nPixels, 1);
distance2 = NaN(nPixels, 1);

distance1(hits==true) = (-b(hits==true) - sqrt(discriminant(hits==true)))./(2*a(hits==true));
distance2(hits==true) = (-b(hits==true) + sqrt(discriminant(hits==true)))./(2*a(hits==true));

distance = min([distance1 distance2],[],2);
distance(b_inside_circle) = max([distance1 distance2],[],2);

hits(distance<0) = false;

distance(hits==false) = NaN;

intersection = zeros(nPixels, 3);
%intersection(hits==true,:) = viewOrigin(hits==true,:) + viewDirection(hits==true,:) .* distance(hits==true);
intersection = viewOrigin + viewDirection .* distance;

%normals = zeros(nPixels, 3);
%normals = normalize(intersection - sphereCenter);
%normals(hits==true,:) = normalize(intersection(hits==true,:) - sphereCenter(hits==true,:));
end