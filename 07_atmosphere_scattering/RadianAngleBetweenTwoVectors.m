function out = RadianAngleBetweenTwoVectors(a,b)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
out = atan2(norm(cross(a,b)), dot(a,b));
end